"""Streaming, coordinate-order-preserving BAM deduplication."""

from __future__ import annotations

import heapq
from collections import deque
from collections.abc import Iterable
from pathlib import Path

import pysam

from .bam_io import atomic_bam_writer, open_coordinate_sorted_bam, validate_paths
from .models import BufferedRecord, Counters, DeduplicationError, DuplicateGroup
from .pairing import PairEvent, ReadEvent, ReleaseEvent, iter_pairing_events


def deduplicate(
    input_path: Path,
    output_path: Path,
    *,
    remove_duplicates: bool = False,
    max_cache_records: int = 1_000_000,
) -> Counters:
    """Deduplicate a coordinate-sorted BAM and atomically publish the output."""
    if max_cache_records < 1:
        raise DeduplicationError("--max-cache-records must be at least 1")
    validate_paths(input_path, output_path)
    counters = Counters()
    with (
        open_coordinate_sorted_bam(input_path) as source,
        atomic_bam_writer(output_path, source.header.to_dict()) as destination,
    ):
        _deduplicate_records(
            source,
            destination,
            counters,
            remove_duplicates=remove_duplicates,
            max_cache_records=max_cache_records,
        )
    return counters


def _deduplicate_records(
    source: Iterable[pysam.AlignedSegment],
    destination: pysam.AlignmentFile,
    counters: Counters,
    *,
    remove_duplicates: bool,
    max_cache_records: int,
) -> None:
    """Manage coordinate ordering and duplicate groups over pairing events."""
    buffer: deque[BufferedRecord] = deque()
    groups: dict[tuple[object, ...], DuplicateGroup] = {}
    group_heap: list[tuple[int, int, tuple[object, ...]]] = []
    group_sequence = 0
    current_reference: int | None = None
    current_start: int | None = None
    last_coordinate: tuple[int, int] | None = None
    seen_unmapped = False

    def flush_ready() -> None:
        while buffer and not buffer[0].blocked:
            item = buffer.popleft()
            if item.drop:
                counters.removed_records += 1
            else:
                destination.write(item.record)

    def resolve_group(group: DuplicateGroup) -> None:
        winner = max(group.pairs, key=lambda pair: (pair.quality, -pair.input_order))
        for pair in group.pairs:
            loser = pair is not winner
            pair.first.record.is_duplicate = loser
            pair.second.record.is_duplicate = loser
            if loser:
                counters.duplicate_pairs += 1
                counters.flagged_records += 2
                if remove_duplicates:
                    pair.first.drop = pair.second.drop = True
            pair.first.blocked = pair.second.blocked = False

    def resolve_before(next_start: int) -> None:
        while group_heap and group_heap[0][0] < next_start:
            _, _, key = heapq.heappop(group_heap)
            group = groups.pop(key, None)
            if group is not None:
                resolve_group(group)

    def resolve_all_groups() -> None:
        while group_heap:
            _, _, key = heapq.heappop(group_heap)
            group = groups.pop(key, None)
            if group is not None:
                resolve_group(group)

    def advance(reference_id: int | None, start: int | None) -> None:
        nonlocal current_reference, current_start
        if current_reference is None:
            current_reference, current_start = reference_id, start
            return
        if reference_id == current_reference and start == current_start:
            return
        if reference_id == current_reference and start is not None:
            resolve_before(start)
        else:
            resolve_all_groups()
        flush_ready()
        current_reference, current_start = reference_id, start

    for event in iter_pairing_events(source):
        if isinstance(event, ReleaseEvent):
            event.item.blocked = False
            counters.passthrough_records += 1
            flush_ready()
            continue
        if isinstance(event, PairEvent):
            pair = event.pair
            counters.eligible_pairs += 1
            group = groups.get(pair.key)
            if group is None:
                group = DuplicateGroup(pair.key, pair.rightmost_start)
                groups[pair.key] = group
                group_sequence += 1
                heapq.heappush(group_heap, (pair.rightmost_start, group_sequence, pair.key))
            group.pairs.append(pair)
            continue

        assert isinstance(event, ReadEvent)
        item = event.item
        record = item.record
        counters.total_records += 1
        if record.is_unmapped:
            seen_unmapped = True
            advance(None, None)
        else:
            coordinate = (record.reference_id, record.reference_start)
            if seen_unmapped or (last_coordinate is not None and coordinate < last_coordinate):
                raise DeduplicationError("input BAM records are not coordinate sorted")
            last_coordinate = coordinate
            advance(*coordinate)

        buffer.append(item)
        counters.peak_cache_records = max(counters.peak_cache_records, len(buffer))
        if len(buffer) > max_cache_records:
            raise DeduplicationError(
                f"cache exceeded --max-cache-records ({max_cache_records}); "
                "increase the limit or use inputs with shorter mate spans"
            )
        if event.passthrough:
            counters.passthrough_records += 1
        flush_ready()

    resolve_all_groups()
    flush_ready()
    if buffer:
        raise DeduplicationError("internal error: records remained blocked at end of input")
