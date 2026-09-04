"""Streaming, coordinate-order-preserving BAM deduplication."""

from __future__ import annotations

import heapq
from collections import deque
from collections.abc import Iterable
from pathlib import Path

import pysam

from .bam_io import atomic_bam_writer, open_coordinate_sorted_bam, validate_paths
from .models import (
    BufferedRecord,
    Counters,
    DeduplicationError,
    DuplicateGroup,
    PairKey,
)
from .pairing import PairEvent, ReadEvent, ReleaseEvent, iter_pairing_events


def deduplicate(
    input_path: Path,
    output_path: Path,
    *,
    remove_duplicates: bool = False,
    max_cache_records: int = 1_000_000,
) -> Counters:
    """Deduplicate a coordinate-sorted BAM and atomically publish the result.

    Lower-quality pairs in an exact ``PairKey`` group receive the SAM duplicate
    flag, or are omitted when ``remove_duplicates`` is set. The active ordered
    cache cannot exceed ``max_cache_records``; invalid options and unsafe input
    raise ``DeduplicationError`` without publishing partial output.
    """
    if max_cache_records < 1:
        raise DeduplicationError("--max-cache-records must be at least 1")
    validate_paths(input_path, output_path)
    counters = Counters()
    with (
        open_coordinate_sorted_bam(input_path) as source,
        atomic_bam_writer(output_path, source.header.to_dict()) as destination,
    ):
        processor = _DeduplicationProcessor(
            destination,
            counters,
            remove_duplicates=remove_duplicates,
            max_cache_records=max_cache_records,
        )
        processor.process(source)
    return counters


class _DeduplicationProcessor:
    """Own the mutable state required for one-pass ordered deduplication."""

    def __init__(
        self,
        destination: pysam.AlignmentFile,
        counters: Counters,
        *,
        remove_duplicates: bool,
        max_cache_records: int,
    ) -> None:
        """Initialize an empty processor with its output and processing options."""
        self.destination = destination
        self.counters = counters
        self.remove_duplicates = remove_duplicates
        self.max_cache_records = max_cache_records
        self.buffer: deque[BufferedRecord] = deque()
        self.groups: dict[PairKey, DuplicateGroup] = {}
        self.group_heap: list[tuple[int, int, PairKey]] = []
        self.group_sequence = 0
        self.current_reference: int | None = None
        self.current_start: int | None = None
        self.last_coordinate: tuple[int, int] | None = None
        self.seen_unmapped = False

    def process(self, source: Iterable[pysam.AlignedSegment]) -> None:
        """Consume pairing events, finalize remaining state, and drain the buffer."""
        for event in iter_pairing_events(source):
            if isinstance(event, ReleaseEvent):
                self._handle_release(event)
            elif isinstance(event, PairEvent):
                self._handle_pair(event)
            else:
                assert isinstance(event, ReadEvent)
                self._handle_read(event)

        self._finalize_groups()
        self._flush_unblocked_records()
        if self.buffer:
            raise DeduplicationError("internal error: records remained blocked at end of input")

    def _handle_release(self, event: ReleaseEvent) -> None:
        """Unblock a non-pairable record and flush any newly writable prefix."""
        event.item.blocked = False
        self.counters.passthrough_records += 1
        self._flush_unblocked_records()

    def _handle_pair(self, event: PairEvent) -> None:
        """Add a completed eligible pair to its duplicate group and expiry heap."""
        pair = event.pair
        self.counters.eligible_pairs += 1
        group = self.groups.get(pair.key)
        if group is None:
            group = DuplicateGroup(pair.key, pair.rightmost_start)
            self.groups[pair.key] = group
            self.group_sequence += 1
            heapq.heappush(
                self.group_heap,
                (pair.rightmost_start, self.group_sequence, pair.key),
            )
        group.pairs.append(pair)

    def _handle_read(self, event: ReadEvent) -> None:
        """Validate, buffer, and account for one input record before it is resolved."""
        item = event.item
        record = item.record
        self.counters.total_records += 1
        self._validate_and_advance_coordinate(record)

        self.buffer.append(item)
        self.counters.peak_cache_records = max(
            self.counters.peak_cache_records, len(self.buffer)
        )
        if len(self.buffer) > self.max_cache_records:
            raise DeduplicationError(
                f"cache exceeded --max-cache-records ({self.max_cache_records}); "
                "increase the limit or use inputs with shorter mate spans"
            )
        if event.passthrough:
            self.counters.passthrough_records += 1
        self._flush_unblocked_records()

    def _validate_and_advance_coordinate(self, record: pysam.AlignedSegment) -> None:
        """Reject unsorted input and move the group-resolution frontier forward."""
        if record.is_unmapped:
            self.seen_unmapped = True
            self._advance_coordinate_frontier(None, None)
            return

        coordinate = (record.reference_id, record.reference_start)
        if self.seen_unmapped or (
            self.last_coordinate is not None and coordinate < self.last_coordinate
        ):
            raise DeduplicationError("input BAM records are not coordinate sorted")
        self.last_coordinate = coordinate
        self._advance_coordinate_frontier(*coordinate)

    def _advance_coordinate_frontier(
        self, reference_id: int | None, start: int | None
    ) -> None:
        """Finalize safe groups and flush after the coordinate stream advances."""
        if self.current_reference is None:
            self.current_reference, self.current_start = reference_id, start
            return
        if reference_id == self.current_reference and start == self.current_start:
            return
        if reference_id == self.current_reference and start is not None:
            self._finalize_groups(start)
        else:
            self._finalize_groups()
        self._flush_unblocked_records()
        self.current_reference, self.current_start = reference_id, start

    def _finalize_groups(self, before_start: int | None = None) -> None:
        """Finalize groups before a cutoff, or all remaining groups when omitted."""
        while self.group_heap and (
            before_start is None or self.group_heap[0][0] < before_start
        ):
            _, _, key = heapq.heappop(self.group_heap)
            group = self.groups.pop(key, None)
            if group is not None:
                self._finalize_duplicate_group(group)

    def _finalize_duplicate_group(self, group: DuplicateGroup) -> None:
        """Select its winner, mark or drop losers, and unblock all pair records."""
        winner = max(
            group.pairs,
            key=lambda pair: (pair.base_quality_sum, -pair.input_order),
        )
        for pair in group.pairs:
            loser = pair is not winner
            pair.first.record.is_duplicate = loser
            pair.second.record.is_duplicate = loser
            if loser:
                self.counters.duplicate_pairs += 1
                self.counters.flagged_records += 2
                if self.remove_duplicates:
                    pair.first.drop = pair.second.drop = True
            pair.first.blocked = pair.second.blocked = False

    def _flush_unblocked_records(self) -> None:
        """Write the unblocked buffer prefix while retaining coordinate order."""
        while self.buffer and not self.buffer[0].blocked:
            item = self.buffer.popleft()
            if item.drop:
                self.counters.removed_records += 1
            else:
                self.destination.write(item.record)
