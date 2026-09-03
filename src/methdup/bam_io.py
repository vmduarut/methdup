"""Validated BAM input and atomic BAM output helpers."""

from __future__ import annotations

import os
import tempfile
from collections.abc import Iterator
from contextlib import contextmanager
from pathlib import Path

import pysam

from .alignment import with_program_record
from .models import DeduplicationError


def validate_paths(input_path: Path, output_path: Path) -> None:
    """Validate the file-path contract shared by the CLI and Python API."""
    if input_path.suffix.lower() != ".bam" or output_path.suffix.lower() != ".bam":
        raise DeduplicationError("input and output paths must both end in .bam")
    if input_path.resolve() == output_path.resolve():
        raise DeduplicationError("input and output paths must be different")


@contextmanager
def open_coordinate_sorted_bam(path: Path) -> Iterator[pysam.AlignmentFile]:
    """Open a BAM and ensure its header declares coordinate sorting."""
    with pysam.AlignmentFile(str(path), "rb") as source:
        if source.header.to_dict().get("HD", {}).get("SO") != "coordinate":
            raise DeduplicationError("input BAM header must declare coordinate sort order")
        yield source


@contextmanager
def atomic_bam_writer(path: Path, header: dict[str, object]) -> Iterator[pysam.AlignmentFile]:
    """Write a BAM to a temporary sibling and publish it only on success."""
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary_name: str | None = None
    try:
        with tempfile.NamedTemporaryFile(
            prefix=f".{path.name}.", suffix=".tmp.bam", dir=path.parent, delete=False
        ) as temporary:
            temporary_name = temporary.name
        with pysam.AlignmentFile(
            temporary_name, "wb", header=with_program_record(header)
        ) as destination:
            yield destination
        os.replace(temporary_name, path)
        temporary_name = None
    finally:
        if temporary_name is not None:
            Path(temporary_name).unlink(missing_ok=True)
