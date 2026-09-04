"""Validated BAM input and atomic BAM output helpers."""

from __future__ import annotations

import os
import tempfile
from collections.abc import Iterator
from contextlib import contextmanager
from pathlib import Path

import pysam

from .alignment import add_program_record
from .models import DeduplicationError


def validate_paths(input_path: Path, output_path: Path) -> None:
    """Validate BAM suffixes and reject an output path equal to the input path.

    Raises ``DeduplicationError`` before opening either file when the path
    contract shared by the CLI and Python API is violated.
    """
    if input_path.suffix.lower() != ".bam" or output_path.suffix.lower() != ".bam":
        raise DeduplicationError("input and output paths must both end in .bam")
    if input_path.resolve() == output_path.resolve():
        raise DeduplicationError("input and output paths must be different")


@contextmanager
def open_coordinate_sorted_bam(path: Path) -> Iterator[pysam.AlignmentFile]:
    """Yield a readable BAM only when its header declares coordinate sorting.

    The caller owns processing inside the context. The file closes on exit;
    a missing ``HD.SO=coordinate`` declaration raises ``DeduplicationError``.
    """
    with pysam.AlignmentFile(str(path), "rb") as source:
        if source.header.to_dict().get("HD", {}).get("SO") != "coordinate":
            raise DeduplicationError("input BAM header must declare coordinate sort order")
        yield source


@contextmanager
def atomic_bam_writer(path: Path, header: dict[str, object]) -> Iterator[pysam.AlignmentFile]:
    """Yield a temporary BAM writer and atomically publish it on success.

    The temporary file is created beside the requested output, receives a
    ``methdup`` program record in its header, and replaces the output path only
    after the writer closes normally. Failures remove the temporary file.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary_name: str | None = None
    try:
        with tempfile.NamedTemporaryFile(
            prefix=f".{path.name}.", suffix=".tmp.bam", dir=path.parent, delete=False
        ) as temporary:
            temporary_name = temporary.name
        with pysam.AlignmentFile(
            temporary_name, "wb", header=add_program_record(header)
        ) as destination:
            yield destination
        os.replace(temporary_name, path)
        temporary_name = None
    finally:
        if temporary_name is not None:
            Path(temporary_name).unlink(missing_ok=True)
