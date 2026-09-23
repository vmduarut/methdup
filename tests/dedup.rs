//! End-to-end deduplication behavior, ported from the original Python suite.

use std::{
    fs::File,
    path::{Path, PathBuf},
    process::Command,
    rc::Rc,
};

use bstr::BString;
use noodles::{
    bam,
    core::Position,
    sam::{
        self,
        alignment::{
            RecordBuf,
            record::{Flags, MappingQuality, cigar::op::Kind},
            record_buf::{Cigar, Sequence},
        },
        header::record::value::{
            Map,
            map::{self, ReadGroup, ReferenceSequence, header::tag as header_tag},
        },
    },
};

fn reference_sequences() -> Vec<(BString, Map<ReferenceSequence>)> {
    (0..2)
        .map(|i| {
            (
                BString::from(format!("chr{}", i + 1)),
                Map::<ReferenceSequence>::new(std::num::NonZero::new(10_000).unwrap()),
            )
        })
        .collect()
}

fn header() -> sam::Header {
    sam::Header::builder()
        .set_header(
            Map::<map::Header>::builder()
                .insert(header_tag::SORT_ORDER, "coordinate")
                .build()
                .expect("valid header map"),
        )
        .set_reference_sequences(reference_sequences().into_iter().collect())
        .build()
}

fn cigar(operations: &str) -> Cigar {
    let mut ops = Vec::new();
    let mut rest = operations;
    loop {
        if rest.is_empty() {
            break;
        }
        let digits = rest.len() - rest.trim_start_matches(|c: char| c.is_ascii_digit()).len();
        let len: usize = rest[..digits].parse().expect("cigar length");
        rest = &rest[digits..];
        let (kind_char, remaining) = rest.split_at(1);
        rest = remaining;
        let kind = match kind_char {
            "M" => Kind::Match,
            "I" => Kind::Insertion,
            "D" => Kind::Deletion,
            "N" => Kind::Skip,
            "S" => Kind::SoftClip,
            "H" => Kind::HardClip,
            "P" => Kind::Pad,
            "=" => Kind::SequenceMatch,
            "X" => Kind::SequenceMismatch,
            other => panic!("unsupported cigar op {other}"),
        };
        ops.push(noodles::sam::alignment::record::cigar::Op::new(kind, len));
    }
    ops.into_iter().collect()
}

fn read(name: &str, start: usize, read1: bool, quality: u8, operations: &str) -> RecordBuf {
    let (flags, mate_start, template_length) = if read1 {
        (
            Flags::SEGMENTED
                | Flags::PROPERLY_SEGMENTED
                | Flags::MATE_REVERSE_COMPLEMENTED
                | Flags::FIRST_SEGMENT,
            200,
            110,
        )
    } else {
        (
            Flags::SEGMENTED
                | Flags::PROPERLY_SEGMENTED
                | Flags::REVERSE_COMPLEMENTED
                | Flags::LAST_SEGMENT,
            100,
            -110,
        )
    };

    RecordBuf::builder()
        .set_name(name)
        .set_flags(flags)
        .set_reference_sequence_id(0)
        .set_alignment_start(Position::try_from(start).expect("1-based start"))
        .set_mapping_quality(MappingQuality::try_from(60).expect("valid mapping quality"))
        .set_cigar(cigar(operations))
        .set_mate_reference_sequence_id(0)
        .set_mate_alignment_start(Position::try_from(mate_start).expect("1-based mate start"))
        .set_template_length(template_length)
        .set_sequence(Sequence::from(&b"AAAAAAAAAA"[..]))
        .set_quality_scores(std::iter::repeat_n(quality, 10).collect())
        .build()
}

fn read_in_rg(
    name: &str,
    start: usize,
    read1: bool,
    quality: u8,
    operations: &str,
    read_group: &str,
) -> RecordBuf {
    use noodles::sam::alignment::{record::data::field::Tag, record_buf::data::field::Value};
    let mut record = read(name, start, read1, quality, operations);
    record
        .data_mut()
        .insert(Tag::READ_GROUP, Value::from(read_group));
    record
}

fn read_group(id: &str, library: &str) -> (BString, Map<ReadGroup>) {
    let map = Map::<ReadGroup>::builder()
        .insert(map::read_group::tag::LIBRARY, library)
        .build()
        .expect("valid read group map");
    (BString::from(id), map)
}

fn header_with_read_groups(read_groups: Vec<(BString, Map<ReadGroup>)>) -> sam::Header {
    let mut header = header();
    header.read_groups_mut().extend(read_groups);
    header
}

fn write_bam_with_header(path: &Path, header: &sam::Header, records: &[RecordBuf]) {
    use noodles::sam::alignment::io::Write as _;
    let mut writer = bam::io::Writer::new(File::create(path).expect("create input bam"));
    writer.write_header(header).expect("write header");
    for record in records {
        writer
            .write_alignment_record(header, record)
            .expect("write record");
    }
    writer.try_finish().expect("finish bam");
}

fn fetch(path: &Path) -> Vec<RecordBuf> {
    let file = File::open(path).expect("open output bam");
    let mut reader = bam::io::Reader::new(file);
    let header = reader.read_header().expect("read header");
    reader
        .record_bufs(&header)
        .collect::<Result<Vec<_>, _>>()
        .expect("read records")
}

fn is_duplicate(record: &RecordBuf) -> bool {
    record.flags().is_duplicate()
}

fn run_dedup(
    input: &Path,
    output: &Path,
    remove_duplicates: bool,
    max_cache_records: u64,
) -> methdup::Counters {
    methdup::deduplicate(input, output, remove_duplicates, max_cache_records)
        .unwrap_or_else(|e| panic!("deduplicate failed: {e}"))
}

/// Every fixture here is written and read back inside one temporary directory.
struct BamFixture {
    _dir: tempfile::TempDir,
    input: PathBuf,
    output: PathBuf,
}

impl BamFixture {
    fn new(records: &[RecordBuf]) -> Self {
        Self::with_header(&header(), records)
    }

    fn with_header(header: &sam::Header, records: &[RecordBuf]) -> Self {
        let _dir = tempfile::tempdir().expect("create temp dir");
        let dir_path = _dir.path().to_path_buf();
        let input = dir_path.join("input.bam");
        let output = dir_path.join("output.bam");
        write_bam_with_header(&input, header, records);
        Self {
            _dir,
            input,
            output,
        }
    }
}

#[test]
fn marks_lower_quality_duplicate_and_keeps_sort_order() {
    let fixture = BamFixture::new(&[
        read("low", 100, true, 10, "10M"),
        read("high", 100, true, 30, "10M"),
        read("low", 200, false, 10, "10M"),
        read("high", 200, false, 30, "10M"),
    ]);
    let counters = run_dedup(&fixture.input, &fixture.output, false, 1_000_000);
    let records = fetch(&fixture.output);
    let names: Vec<_> = records
        .iter()
        .map(|record| record.name().unwrap().to_string())
        .collect();
    assert_eq!(names, ["low", "high", "low", "high"]);
    let duplicates: Vec<_> = records.iter().map(is_duplicate).collect();
    assert_eq!(duplicates, [true, false, true, false]);
    assert_eq!(
        (
            counters.eligible_pairs,
            counters.duplicate_pairs,
            counters.flagged_records
        ),
        (2, 1, 2)
    );
    let file = File::open(&fixture.output).expect("open output bam");
    let mut reader = bam::io::Reader::new(file);
    let header = reader.read_header().expect("read header");
    let program_names: Vec<_> = header
        .programs()
        .as_ref()
        .keys()
        .map(|id| id.to_string())
        .collect();
    assert_eq!(program_names.last().map(String::as_str), Some("methdup"));
}

#[test]
fn remove_mode_omits_both_losing_records() {
    let fixture = BamFixture::new(&[
        read("low", 100, true, 10, "10M"),
        read("high", 100, true, 30, "10M"),
        read("low", 200, false, 10, "10M"),
        read("high", 200, false, 30, "10M"),
    ]);
    let counters = run_dedup(&fixture.input, &fixture.output, true, 1_000_000);
    let records = fetch(&fixture.output);
    let names: Vec<_> = records
        .iter()
        .map(|record| record.name().unwrap().to_string())
        .collect();
    assert_eq!(names, ["high", "high"]);
    assert_eq!(counters.removed_records, 2);
}

#[test]
fn cigar_changes_duplicate_identity() {
    let fixture = BamFixture::new(&[
        read("a", 100, true, 10, "10M"),
        read("b", 100, true, 30, "10M"),
        read("a", 200, false, 10, "10M"),
        read("b", 200, false, 30, "5M1I4M"),
    ]);
    run_dedup(&fixture.input, &fixture.output, false, 1_000_000);
    let records = fetch(&fixture.output);
    assert!(!records.iter().any(is_duplicate));
}

#[test]
fn recomputes_existing_duplicate_flag_and_uses_input_order_for_ties() {
    let first1 = read("first", 100, true, 20, "10M");
    let first2 = read("first", 200, false, 20, "10M");
    let mut second1 = read("second", 100, true, 20, "10M");
    let mut second2 = read("second", 200, false, 20, "10M");
    second1.flags_mut().insert(Flags::DUPLICATE);
    second2.flags_mut().insert(Flags::DUPLICATE);

    let fixture = BamFixture::new(&[first1, second1, second2, first2]);
    run_dedup(&fixture.input, &fixture.output, false, 1_000_000);
    let records = fetch(&fixture.output);
    let duplicates: Vec<_> = records.iter().map(is_duplicate).collect();
    assert_eq!(duplicates, [false, true, true, false]);
}

#[test]
fn cache_limit_on_active_window_does_not_publish_output() {
    let fixture = BamFixture::new(&[
        read("a", 100, true, 20, "10M"),
        read("pass", 101, true, 20, "10M"),
    ]);
    let result = methdup::deduplicate(&fixture.input, &fixture.output, false, 1);
    let message = result.expect_err("cache limit exceeded").to_string();
    assert!(message.starts_with("cache exceeded"), "{message}");
    assert!(!fixture.output.exists());
}

#[test]
fn incomplete_and_cross_reference_pairs_pass_through_unchanged() {
    let incomplete = read("incomplete", 100, true, 20, "10M");
    let mut cross_reference = read("cross", 101, true, 20, "10M");
    cross_reference.mate_reference_sequence_id_mut().replace(1);

    let fixture = BamFixture::new(&[incomplete, cross_reference]);
    let counters = run_dedup(&fixture.input, &fixture.output, false, 1_000_000);
    let records = fetch(&fixture.output);
    let names: Vec<_> = records
        .iter()
        .map(|record| record.name().unwrap().to_string())
        .collect();
    assert_eq!(names, ["incomplete", "cross"]);
    assert!(!records.iter().any(is_duplicate));
    assert_eq!(counters.passthrough_records, 2);
}

#[test]
fn rejects_unsorted_input_without_publishing_output() {
    let fixture = BamFixture::new(&[
        read("a", 200, false, 20, "10M"),
        read("a", 100, true, 20, "10M"),
    ]);
    let result = methdup::deduplicate(&fixture.input, &fixture.output, false, 1_000_000);
    let message = result.expect_err("unsorted input rejected").to_string();
    assert!(message.contains("not coordinate sorted"), "{message}");
    assert!(!fixture.output.exists());
}

#[test]
fn cache_limit_fails_safely_and_cli_reports_error() {
    let fixture = BamFixture::new(&[read("a", 100, true, 20, "10M")]);
    let output_text = Command::new(env!("CARGO_BIN_EXE_methdup"))
        .arg(&fixture.input)
        .arg(&fixture.output)
        .arg("--max-cache-records")
        .arg("0")
        .output()
        .expect("run cli");
    assert_eq!(output_text.status.code(), Some(2));
    let stderr = String::from_utf8_lossy(&output_text.stderr);
    assert!(stderr.contains("must be at least 1"), "{stderr}");
    assert!(!fixture.output.exists());
}

#[test]
fn cli_reports_summary_on_success() {
    let fixture = BamFixture::new(&[
        read("low", 100, true, 10, "10M"),
        read("high", 100, true, 30, "10M"),
        read("low", 200, false, 10, "10M"),
        read("high", 200, false, 30, "10M"),
    ]);
    let output_text = Command::new(env!("CARGO_BIN_EXE_methdup"))
        .arg(&fixture.input)
        .arg(&fixture.output)
        .output()
        .expect("run cli");
    assert_eq!(
        output_text.status.code(),
        Some(0),
        "{}",
        String::from_utf8_lossy(&output_text.stderr)
    );
    let stderr = String::from_utf8_lossy(&output_text.stderr);
    assert!(
        stderr.contains("records=4 eligible_pairs=2 duplicate_pairs=1 flagged_records=2"),
        "{stderr}"
    );
    assert!(fixture.output.exists());
}

#[test]
fn pairing_emits_reads_before_completed_pair() {
    use methdup::pairing::{Pairing, PairingEvent};

    let header = header();
    let mut pairing = Pairing::new(&header);
    let mut events = pairing.push(read("pair", 100, true, 20, "10M"));
    assert_eq!(events.len(), 1);
    let first = match events.pop() {
        Some(PairingEvent::Read { item, passthrough }) => {
            assert!(!passthrough);
            item
        }
        other => panic!("expected ReadEvent, got {other:?}"),
    };

    let events = pairing.push(read("pair", 200, false, 20, "10M"));
    assert_eq!(events.len(), 2);
    let second = match &events[0] {
        PairingEvent::Read { item, .. } => item.clone(),
        other => panic!("expected ReadEvent, got {other:?}"),
    };
    match &events[1] {
        PairingEvent::Pair {
            first: pair_first,
            second: pair_second,
            ..
        } => {
            assert!(Rc::ptr_eq(pair_first, &first));
            assert!(Rc::ptr_eq(pair_second, &second));
        }
        other => panic!("expected PairEvent, got {other:?}"),
    }
}

#[test]
fn pairing_emits_release_for_incomplete_pair() {
    use methdup::pairing::{Pairing, PairingEvent};

    let header = header();
    let mut pairing = Pairing::new(&header);
    let events = pairing.push(read("incomplete", 100, true, 20, "10M"));
    assert_eq!(events.len(), 1);
    let item = match &events[0] {
        PairingEvent::Read { item, .. } => item.clone(),
        other => panic!("expected ReadEvent, got {other:?}"),
    };
    let events = pairing.finish();
    assert_eq!(events.len(), 1);
    match &events[0] {
        PairingEvent::Release { item: released } => {
            assert!(Rc::ptr_eq(released, &item));
        }
        other => panic!("expected ReleaseEvent, got {other:?}"),
    }
}

#[test]
fn different_libraries_are_not_duplicates() {
    let header =
        header_with_read_groups(vec![read_group("lib1", "lib1"), read_group("lib2", "lib2")]);
    let fixture = BamFixture::with_header(
        &header,
        &[
            read_in_rg("a", 100, true, 10, "10M", "lib1"),
            read_in_rg("b", 100, true, 30, "10M", "lib2"),
            read_in_rg("a", 200, false, 10, "10M", "lib1"),
            read_in_rg("b", 200, false, 30, "10M", "lib2"),
        ],
    );
    let counters = run_dedup(&fixture.input, &fixture.output, false, 1_000_000);
    let records = fetch(&fixture.output);
    assert!(!records.iter().any(is_duplicate));
    assert_eq!((counters.eligible_pairs, counters.duplicate_pairs), (2, 0));
}

#[test]
fn different_read_groups_share_duplicate_identity() {
    let header =
        header_with_read_groups(vec![read_group("rg1", "libX"), read_group("rg2", "libX")]);
    let fixture = BamFixture::with_header(
        &header,
        &[
            read_in_rg("low", 100, true, 10, "10M", "rg1"),
            read_in_rg("high", 100, true, 30, "10M", "rg2"),
            read_in_rg("low", 200, false, 10, "10M", "rg1"),
            read_in_rg("high", 200, false, 30, "10M", "rg2"),
        ],
    );
    let counters = run_dedup(&fixture.input, &fixture.output, false, 1_000_000);
    let records = fetch(&fixture.output);
    let duplicates: Vec<_> = records.iter().map(is_duplicate).collect();
    assert_eq!(duplicates, [true, false, true, false]);
    assert_eq!((counters.eligible_pairs, counters.duplicate_pairs), (2, 1));
}

#[test]
fn missing_read_group_is_distinct_from_named_library() {
    let header = header_with_read_groups(vec![read_group("rg1", "lib1")]);
    let fixture = BamFixture::with_header(
        &header,
        &[
            read_in_rg("a", 100, true, 10, "10M", "rg1"),
            read("b", 100, true, 30, "10M"),
            read_in_rg("a", 200, false, 10, "10M", "rg1"),
            read("b", 200, false, 30, "10M"),
        ],
    );
    run_dedup(&fixture.input, &fixture.output, false, 1_000_000);
    let records = fetch(&fixture.output);
    assert!(!records.iter().any(is_duplicate));
}

#[test]
fn same_name_in_different_libraries_do_not_cross_pair() {
    let header =
        header_with_read_groups(vec![read_group("lib1", "lib1"), read_group("lib2", "lib2")]);
    let fixture = BamFixture::with_header(
        &header,
        &[
            read_in_rg("x", 100, true, 10, "10M", "lib1"),
            read_in_rg("x", 100, true, 30, "10M", "lib2"),
            read_in_rg("x", 200, false, 20, "10M", "lib1"),
            read_in_rg("x", 200, false, 40, "10M", "lib2"),
        ],
    );
    let counters = run_dedup(&fixture.input, &fixture.output, false, 1_000_000);
    let records = fetch(&fixture.output);
    let duplicates: Vec<_> = records.iter().map(is_duplicate).collect();
    assert_eq!(duplicates, [false, false, false, false]);
    assert_eq!((counters.eligible_pairs, counters.duplicate_pairs), (2, 0));
}
