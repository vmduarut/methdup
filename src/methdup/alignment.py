"""BAM alignment predicates and duplicate-key construction."""

from __future__ import annotations

import pysam

from .models import PairKey


def is_candidate(record: pysam.AlignedSegment) -> bool:
    """Return whether a record can be one end of an eligible duplicate pair."""
    return (
        record.is_paired
        and record.is_proper_pair
        and not record.is_unmapped
        and not record.mate_is_unmapped
        and not record.is_secondary
        and not record.is_supplementary
        and (record.is_read1 != record.is_read2)
        and record.reference_id == record.next_reference_id
        and record.next_reference_start >= 0
    )


def pair_key(read1: pysam.AlignedSegment, read2: pysam.AlignedSegment) -> PairKey:
    """Return the full alignment identity used to group duplicate fragments."""
    return PairKey(
        read1.reference_id,
        read1.reference_start,
        read1.is_reverse,
        read1.cigarstring,
        read2.reference_id,
        read2.reference_start,
        read2.is_reverse,
        read2.cigarstring,
    )


def are_reciprocal(read1: pysam.AlignedSegment, read2: pysam.AlignedSegment) -> bool:
    """Return whether each mate's declared coordinates point to the other."""
    return (
        read1.next_reference_id == read2.reference_id
        and read2.next_reference_id == read1.reference_id
        and read1.next_reference_start == read2.reference_start
        and read2.next_reference_start == read1.reference_start
    )


def with_program_record(header: dict[str, object]) -> dict[str, object]:
    """Copy a BAM header and append a unique methdup program record."""
    result = dict(header)
    programs = list(result.get("PG", []))
    existing_ids = {entry.get("ID") for entry in programs if isinstance(entry, dict)}
    program_id, suffix = "methdup", 1
    while program_id in existing_ids:
        suffix += 1
        program_id = f"methdup.{suffix}"
    programs.append({"ID": program_id, "PN": "methdup", "VN": "0.1.0", "CL": "methdup"})
    result["PG"] = programs
    return result
