from __future__ import annotations

from pathlib import Path

import pysam
import pytest

from methdup import DeduplicationError, deduplicate, run
from methdup.pairing import PairEvent, ReadEvent, ReleaseEvent, iter_pairing_events

HEADER = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": "chr1", "LN": 10_000}]}


def read(name: str, *, start: int, read1: bool, quality: int, cigar: str = "10M") -> pysam.AlignedSegment:
    record = pysam.AlignedSegment()
    record.query_name = name
    record.query_sequence = "A" * 10
    record.flag = 99 if read1 else 147
    record.reference_id = 0
    record.reference_start = start
    record.mapping_quality = 60
    record.cigarstring = cigar
    record.next_reference_id = 0
    record.next_reference_start = 200 if read1 else 100
    record.template_length = 110 if read1 else -110
    record.query_qualities = pysam.qualitystring_to_array(chr(quality + 33) * 10)
    return record


def write_bam(path: Path, records: list[pysam.AlignedSegment], header: dict = HEADER) -> None:
    with pysam.AlignmentFile(path, "wb", header=header) as bam:
        for record in records:
            bam.write(record)


def fetch(path: Path) -> list[pysam.AlignedSegment]:
    with pysam.AlignmentFile(path, "rb") as bam:
        return list(bam)


def test_marks_lower_quality_duplicate_and_keeps_sort_order(tmp_path: Path) -> None:
    input_bam, output_bam = tmp_path / "input.bam", tmp_path / "output.bam"
    write_bam(input_bam, [read("low", start=100, read1=True, quality=10), read("high", start=100, read1=True, quality=30), read("low", start=200, read1=False, quality=10), read("high", start=200, read1=False, quality=30)])
    counters = deduplicate(input_bam, output_bam)
    records = fetch(output_bam)
    assert [record.query_name for record in records] == ["low", "high", "low", "high"]
    assert [record.is_duplicate for record in records] == [True, False, True, False]
    assert (counters.eligible_pairs, counters.duplicate_pairs, counters.flagged_records) == (2, 1, 2)
    with pysam.AlignmentFile(output_bam, "rb") as bam:
        assert bam.header.to_dict()["PG"][-1]["PN"] == "methdup"


def test_remove_mode_omits_both_losing_records(tmp_path: Path) -> None:
    input_bam, output_bam = tmp_path / "input.bam", tmp_path / "output.bam"
    write_bam(input_bam, [read("low", start=100, read1=True, quality=10), read("high", start=100, read1=True, quality=30), read("low", start=200, read1=False, quality=10), read("high", start=200, read1=False, quality=30)])
    counters = deduplicate(input_bam, output_bam, remove_duplicates=True)
    assert [record.query_name for record in fetch(output_bam)] == ["high", "high"]
    assert counters.removed_records == 2


def test_cigar_changes_duplicate_identity(tmp_path: Path) -> None:
    input_bam, output_bam = tmp_path / "input.bam", tmp_path / "output.bam"
    write_bam(input_bam, [read("a", start=100, read1=True, quality=10), read("b", start=100, read1=True, quality=30), read("a", start=200, read1=False, quality=10), read("b", start=200, read1=False, quality=30, cigar="5M1I4M")])
    deduplicate(input_bam, output_bam)
    assert not any(record.is_duplicate for record in fetch(output_bam))


def test_recomputes_existing_duplicate_flag_and_uses_input_order_for_ties(tmp_path: Path) -> None:
    input_bam, output_bam = tmp_path / "input.bam", tmp_path / "output.bam"
    first1, first2 = read("first", start=100, read1=True, quality=20), read("first", start=200, read1=False, quality=20)
    second1, second2 = read("second", start=100, read1=True, quality=20), read("second", start=200, read1=False, quality=20)
    second1.is_duplicate = second2.is_duplicate = True
    write_bam(input_bam, [first1, second1, second2, first2])
    deduplicate(input_bam, output_bam)
    records = fetch(output_bam)
    assert [record.is_duplicate for record in records] == [False, True, True, False]


def test_cache_limit_on_active_window_does_not_publish_output(tmp_path: Path) -> None:
    input_bam, output_bam = tmp_path / "input.bam", tmp_path / "output.bam"
    write_bam(input_bam, [read("a", start=100, read1=True, quality=20), read("pass", start=101, read1=True, quality=20)])
    with pytest.raises(DeduplicationError, match="cache exceeded"):
        deduplicate(input_bam, output_bam, max_cache_records=1)
    assert not output_bam.exists()


def test_incomplete_and_cross_reference_pairs_pass_through_unchanged(tmp_path: Path) -> None:
    input_bam, output_bam = tmp_path / "input.bam", tmp_path / "output.bam"
    incomplete = read("incomplete", start=100, read1=True, quality=20)
    cross_reference = read("cross", start=101, read1=True, quality=20)
    cross_reference.next_reference_id = 1
    header = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": "chr1", "LN": 10_000}, {"SN": "chr2", "LN": 10_000}]}
    write_bam(input_bam, [incomplete, cross_reference], header)
    counters = deduplicate(input_bam, output_bam)
    records = fetch(output_bam)
    assert [record.query_name for record in records] == ["incomplete", "cross"]
    assert not any(record.is_duplicate for record in records)
    assert counters.passthrough_records == 2


def test_rejects_unsorted_input_without_publishing_output(tmp_path: Path) -> None:
    input_bam, output_bam = tmp_path / "input.bam", tmp_path / "output.bam"
    write_bam(input_bam, [read("a", start=200, read1=False, quality=20), read("a", start=100, read1=True, quality=20)])
    with pytest.raises(DeduplicationError, match="not coordinate sorted"):
        deduplicate(input_bam, output_bam)
    assert not output_bam.exists()


def test_cache_limit_fails_safely_and_cli_reports_error(tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
    input_bam, output_bam = tmp_path / "input.bam", tmp_path / "output.bam"
    write_bam(input_bam, [read("a", start=100, read1=True, quality=20)])
    assert run([str(input_bam), str(output_bam), "--max-cache-records", "0"]) == 2
    assert "must be at least 1" in capsys.readouterr().err
    assert not output_bam.exists()


def test_pairing_generator_emits_reads_before_completed_pair() -> None:
    events = list(
        iter_pairing_events(
            [read("pair", start=100, read1=True, quality=20), read("pair", start=200, read1=False, quality=20)]
        )
    )
    assert [type(event) for event in events] == [ReadEvent, ReadEvent, PairEvent]
    assert events[0].item.record.is_read1
    assert events[1].item.record.is_read2
    assert events[2].pair.first is events[0].item
    assert events[2].pair.second is events[1].item


def test_pairing_generator_emits_release_for_incomplete_pair() -> None:
    events = list(iter_pairing_events([read("incomplete", start=100, read1=True, quality=20)]))
    assert [type(event) for event in events] == [ReadEvent, ReleaseEvent]
    assert events[1].item is events[0].item


def test_pairing_generator_marks_unmapped_record_as_immediate_passthrough() -> None:
    unmapped = read("unmapped", start=100, read1=True, quality=20)
    unmapped.flag = 4
    events = list(iter_pairing_events([unmapped]))
    assert [type(event) for event in events] == [ReadEvent]
    assert events[0].passthrough is True
