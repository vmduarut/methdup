//! BAM alignment predicates and duplicate-key construction.

use bstr::BString;
use noodles::sam::{
    Header,
    alignment::{
        RecordBuf,
        record::cigar::{Cigar as _, op::Kind},
        record::{self, Flags},
        record_buf::{Cigar, data::field::Value},
    },
    header::record::value::{
        Map,
        map::{self, program::tag, read_group::tag::LIBRARY},
    },
};

use crate::models::{DeduplicationError, PairKey};

/// Fallback library name when a record has no `RG` tag, its read group is
/// missing from the header, or the read group declares no library.
const UNKNOWN_LIBRARY: &str = "Unknown Library";

/// Returns whether an alignment is eligible for duplicate-pair matching.
///
/// Eligible records are mapped, primary, same-reference proper-pair ends with
/// exactly one read-end flag and a known mate coordinate. All other records
/// pass through the deduplication pipeline unchanged.
pub fn is_eligible_pair_end(record: &RecordBuf) -> bool {
    let flags = record.flags();
    let same_reference = record
        .reference_sequence_id()
        .is_some_and(|r| record.mate_reference_sequence_id() == Some(r));
    flags.is_segmented()
        && flags.is_properly_segmented()
        && !flags.is_unmapped()
        && !flags.is_mate_unmapped()
        && !flags.is_secondary()
        && !flags.is_supplementary()
        && flags.is_first_segment() != flags.is_last_segment()
        && same_reference
        && record.mate_alignment_start().is_some()
}

/// Converts a record's CIGAR into its SAM string form, if non-empty.
fn cigar_string(cigar: &Cigar) -> Option<String> {
    if cigar.is_empty() {
        return None;
    }
    Some(cigar_to_string(cigar))
}

fn op_kind_char(kind: Kind) -> char {
    match kind {
        Kind::Match => 'M',
        Kind::Insertion => 'I',
        Kind::Deletion => 'D',
        Kind::Skip => 'N',
        Kind::SoftClip => 'S',
        Kind::HardClip => 'H',
        Kind::Pad => 'P',
        Kind::SequenceMatch => '=',
        Kind::SequenceMismatch => 'X',
    }
}

/// Formats a CIGAR as a SAM cigar string (e.g. `10M2I8M`).
fn cigar_to_string(cigar: &Cigar) -> String {
    let mut string = String::new();
    for op in cigar.as_ref() {
        string.push_str(&op.len().to_string());
        string.push(op_kind_char(op.kind()));
    }
    string
}

/// Resolves a record's library by following its `RG` tag into the header's
/// `@RG` records and reading the `LB` field.
///
/// Records with no `RG`, an unknown read group, or no `LB` fall back to a
/// single shared unknown-library group so they are compared with each other.
pub(crate) fn record_library(header: &Header, record: &RecordBuf) -> BString {
    if let Some(Value::String(read_group_id)) =
        record.data().get(&record::data::field::Tag::READ_GROUP)
        && let Some(read_group) = header.read_groups().get(read_group_id)
        && let Some(library) = read_group.other_fields().get(&LIBRARY)
    {
        return library.clone();
    }
    BString::from(UNKNOWN_LIBRARY)
}

/// Builds the exact alignment identity used to group duplicate fragments.
pub fn build_pair_key(read1: &RecordBuf, read2: &RecordBuf, header: &Header) -> PairKey {
    PairKey {
        read1_reference: reference_id(read1),
        read1_start: reference_start(read1),
        read1_reverse: read1.flags().is_reverse_complemented(),
        read1_cigar: cigar_string(read1.cigar()),
        read2_reference: reference_id(read2),
        read2_start: reference_start(read2),
        read2_reverse: read2.flags().is_reverse_complemented(),
        read2_cigar: cigar_string(read2.cigar()),
        library: record_library(header, read1),
    }
}

/// Returns whether two records' mate coordinates point to one another.
///
/// This guards against assembling a same-name read1/read2 pair whose SAM mate
/// fields describe different alignments.
pub fn have_reciprocal_mate_coordinates(read1: &RecordBuf, read2: &RecordBuf) -> bool {
    read1.mate_reference_sequence_id() == read2.reference_sequence_id()
        && read2.mate_reference_sequence_id() == read1.reference_sequence_id()
        && mate_start(read1) == reference_start(read2)
        && mate_start(read2) == reference_start(read1)
}

/// Returns a mapped record's reference sequence ID as a signed value.
pub fn reference_id(record: &RecordBuf) -> i32 {
    record.reference_sequence_id().map_or(-1, |r| r as i32)
}

/// Returns a mapped record's zero-based reference start as a signed value.
pub fn reference_start(record: &RecordBuf) -> i32 {
    record
        .alignment_start()
        .map_or(-1, |p| i32::try_from(p.get() - 1).unwrap_or(i32::MAX))
}

/// Returns a record's zero-based mate reference start, or -1 when unset.
pub fn mate_start(record: &RecordBuf) -> i32 {
    record
        .mate_alignment_start()
        .map_or(-1, |p| i32::try_from(p.get() - 1).unwrap_or(i32::MAX))
}

/// Adds a unique `methdup` program record to the BAM header.
///
/// Existing program records are preserved. If `methdup` is already an ID, a
/// numeric suffix is added so the output header remains valid.
pub fn add_program_record(header: &mut noodles::sam::Header) -> Result<(), DeduplicationError> {
    const PROGRAM_NAME: &str = "methdup";
    const PROGRAM_VERSION: &str = env!("CARGO_PKG_VERSION");

    let existing_ids = header
        .programs()
        .as_ref()
        .keys()
        .cloned()
        .collect::<Vec<_>>();

    let mut program_id = BString::from(PROGRAM_NAME);
    let mut suffix = 1usize;
    while existing_ids.contains(&program_id) {
        suffix += 1;
        program_id = BString::from(format!("{PROGRAM_NAME}.{suffix}"));
    }

    let program = Map::<map::Program>::builder()
        .insert(tag::NAME, PROGRAM_NAME)
        .insert(tag::VERSION, PROGRAM_VERSION)
        .insert(tag::COMMAND_LINE, PROGRAM_NAME)
        .build()
        .map_err(|e| DeduplicationError::Other(e.to_string()))?;

    header.programs_mut().as_mut().insert(program_id, program);

    Ok(())
}

/// Sets or clears the SAM duplicate flag (`0x400`) on a record.
pub fn set_duplicate(record: &mut RecordBuf, is_duplicate: bool) {
    if is_duplicate {
        record.flags_mut().insert(Flags::DUPLICATE);
    } else {
        record.flags_mut().remove(Flags::DUPLICATE);
    }
}
