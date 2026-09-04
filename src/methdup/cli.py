"""Command-line parsing and terminal reporting for methdup."""

from __future__ import annotations

import argparse
import sys
from collections.abc import Sequence
from pathlib import Path
from typing import TextIO

from .models import Counters, DeduplicationError
from .processor import deduplicate


def build_parser() -> argparse.ArgumentParser:
    """Build the command-line parser for BAM paths and deduplication options."""
    parser = argparse.ArgumentParser(description="One-pass duplicate marking for coordinate-sorted BAM files.")
    parser.add_argument("input_bam", type=Path)
    parser.add_argument("output_bam", type=Path)
    parser.add_argument("--remove-duplicates", action="store_true", help="omit losing pairs instead of writing them with SAM flag 0x400")
    parser.add_argument("--max-cache-records", type=int, default=1_000_000, help="maximum buffered records before failing safely")
    return parser


def format_processing_summary(counters: Counters) -> str:
    """Format all processing counters as the single success-line CLI report."""
    return (
        "methdup: "
        f"records={counters.total_records} eligible_pairs={counters.eligible_pairs} "
        f"duplicate_pairs={counters.duplicate_pairs} flagged_records={counters.flagged_records} "
        f"removed_records={counters.removed_records} passthrough_records={counters.passthrough_records} "
        f"peak_cache_records={counters.peak_cache_records}"
    )


def run(argv: Sequence[str] | None = None, *, stderr: TextIO | None = None) -> int:
    """Run the CLI, report its result to stderr, and return a process status.

    Argument parsing errors are handled by ``argparse``. Processing and file
    errors produce a prefixed diagnostic and status 2; success prints counters
    and returns 0. ``stderr`` supports injection by callers and tests.
    """
    args = build_parser().parse_args(argv)
    stderr = sys.stderr if stderr is None else stderr
    try:
        counters = deduplicate(args.input_bam, args.output_bam, remove_duplicates=args.remove_duplicates, max_cache_records=args.max_cache_records)
    except (DeduplicationError, OSError, ValueError) as error:
        print(f"methdup: error: {error}", file=stderr)
        return 2
    print(format_processing_summary(counters), file=stderr)
    return 0


def main() -> None:
    """Exit the process using the status code returned by :func:`run`."""
    raise SystemExit(run())
