"""BAM alignment predicates and duplicate-key construction."""

from __future__ import annotations

import pysam

from .models import PairKey


def is_eligible_pair_end(record: pysam.AlignedSegment) -> bool:
    """Return whether an alignment is eligible for duplicate-pair matching.

    Eligible records are mapped, primary, same-reference proper-pair ends with
    exactly one read-end flag and a known mate coordinate. All other records
    pass through the deduplication pipeline unchanged.
    """
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


def build_pair_key(read1: pysam.AlignedSegment, read2: pysam.AlignedSegment) -> PairKey:
    """Build the exact alignment identity used to group duplicate fragments.

    The key keeps both ends in read1/read2 order and includes each end's
    reference, start, strand, and CIGAR so alignments with different internal
    structures are not collapsed into the same duplicate group.
    """
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


def have_reciprocal_mate_coordinates(
    read1: pysam.AlignedSegment, read2: pysam.AlignedSegment
) -> bool:
    """Return whether two records' mate coordinates point to one another.

    This guards against assembling a same-name read1/read2 pair whose SAM mate
    fields describe different alignments.
    """
    return (
        read1.next_reference_id == read2.reference_id
        and read2.next_reference_id == read1.reference_id
        and read1.next_reference_start == read2.reference_start
        and read2.next_reference_start == read1.reference_start
    )


def add_program_record(header: dict[str, object]) -> dict[str, object]:
    """Return a copied BAM header with a unique ``methdup`` program entry.

    Existing program records are preserved. If ``methdup`` is already an ID,
    a numeric suffix is added so the output header remains valid.
    """
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
