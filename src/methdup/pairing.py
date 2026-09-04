"""Streaming assembly of primary paired-end alignments."""

from __future__ import annotations

from collections.abc import Iterable, Iterator
from dataclasses import dataclass

import pysam

from .alignment import (
    build_pair_key,
    have_reciprocal_mate_coordinates,
    is_eligible_pair_end,
)
from .models import BufferedRecord, Pair


@dataclass
class ReadEvent:
    """An input alignment that must be added to the ordered output buffer."""

    item: BufferedRecord
    passthrough: bool = False


@dataclass
class PairEvent:
    """A complete reciprocal pair whose two reads were already emitted."""

    pair: Pair


@dataclass
class ReleaseEvent:
    """An eligible record which cannot form a valid pair and may be written."""

    item: BufferedRecord


type PairingEvent = ReadEvent | PairEvent | ReleaseEvent


def iter_pairing_events(records: Iterable[pysam.AlignedSegment]) -> Iterator[PairingEvent]:
    """Stream records into read, completed-pair, and safe-release events.

    Every input record first produces a ``ReadEvent`` for ordered buffering.
    Eligible mates are retained by query name until a reciprocal read1/read2
    pair completes, or released when coordinate order proves a mate cannot
    arrive. Invalid, incomplete, and non-candidate records are released for
    unchanged output.
    """
    pending: dict[str, dict[int, BufferedRecord]] = {}
    invalid_names: set[str] = set()
    record_sequence = 0
    current_reference: int | None = None
    current_start: int | None = None

    def release_expired_pending(
        next_reference: int | None, next_start: int | None
    ) -> list[BufferedRecord]:
        """Remove pending records whose declared mate position is now behind us.

        A reference transition or unmapped input releases every pending record;
        otherwise a record is released once the coordinate stream has passed
        its mate's expected start. Returned objects are the same buffered
        records already emitted in earlier ``ReadEvent`` instances.
        """
        released: list[BufferedRecord] = []
        for query_name, ends in list(pending.items()):
            for end, item in list(ends.items()):
                record = item.record
                passed_mate = (
                    next_reference is None
                    or record.next_reference_id != next_reference
                    or (next_start is not None and record.next_reference_start < next_start)
                )
                if passed_mate:
                    released.append(item)
                    del ends[end]
            if not ends:
                del pending[query_name]
        return released

    for record in records:
        if record.is_unmapped:
            released = release_expired_pending(None, None)
            invalid_names = set()
            current_reference, current_start = None, None
        elif current_reference != record.reference_id:
            released = release_expired_pending(None, None)
            invalid_names = set()
            current_reference, current_start = record.reference_id, record.reference_start
        elif current_start != record.reference_start:
            released = release_expired_pending(record.reference_id, record.reference_start)
            current_start = record.reference_start
        else:
            released = []
        for item in released:
            yield ReleaseEvent(item)

        record_sequence += 1
        candidate = is_eligible_pair_end(record)
        item = BufferedRecord(record, blocked=candidate, ordinal=record_sequence)
        yield ReadEvent(item, passthrough=not candidate)
        if not candidate:
            continue
        if record.query_name in invalid_names:
            yield ReleaseEvent(item)
            continue

        end, other_end = (1, 2) if record.is_read1 else (2, 1)
        ends = pending.setdefault(record.query_name, {})
        if end in ends:
            previous = ends.pop(end)
            if not ends:
                del pending[record.query_name]
            invalid_names.add(record.query_name)
            yield ReleaseEvent(previous)
            yield ReleaseEvent(item)
            continue

        other = ends.get(other_end)
        if other is None:
            if current_start is not None and record.next_reference_start < current_start:
                if not ends:
                    del pending[record.query_name]
                yield ReleaseEvent(item)
                continue
            ends[end] = item
            continue

        del ends[other_end]
        if not ends:
            del pending[record.query_name]
        read1_item, read2_item = (item, other) if end == 1 else (other, item)
        if not have_reciprocal_mate_coordinates(read1_item.record, read2_item.record):
            yield ReleaseEvent(read1_item)
            yield ReleaseEvent(read2_item)
            continue
        yield PairEvent(
            Pair(
                read1_item,
                read2_item,
                build_pair_key(read1_item.record, read2_item.record),
                max(read1_item.record.reference_start, read2_item.record.reference_start),
                min(read1_item.ordinal, read2_item.ordinal),
            )
        )

    for item in release_expired_pending(None, None):
        yield ReleaseEvent(item)
