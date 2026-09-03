"""Data structures shared by the duplicate-detection pipeline."""

from __future__ import annotations

from dataclasses import dataclass, field

import pysam


class DeduplicationError(RuntimeError):
    """Raised when input cannot safely be deduplicated in one pass."""


@dataclass
class Counters:
    total_records: int = 0
    eligible_pairs: int = 0
    duplicate_pairs: int = 0
    flagged_records: int = 0
    removed_records: int = 0
    passthrough_records: int = 0
    peak_cache_records: int = 0


@dataclass
class BufferedRecord:
    record: pysam.AlignedSegment
    blocked: bool
    ordinal: int
    drop: bool = False


@dataclass
class Pair:
    first: BufferedRecord
    second: BufferedRecord
    key: tuple[object, ...]
    rightmost_start: int
    input_order: int

    @property
    def quality(self) -> int:
        return sum(self.first.record.query_qualities or ()) + sum(self.second.record.query_qualities or ())


@dataclass
class DuplicateGroup:
    key: tuple[object, ...]
    rightmost_start: int
    pairs: list[Pair] = field(default_factory=list)
