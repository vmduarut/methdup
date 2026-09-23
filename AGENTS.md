# Repository Instructions

## Commands

- Use the latest stable Rust; `rust-toolchain.toml` selects the stable channel and required components for rustup-based environments. Run `rustup update stable` when rustup is available. The crate uses edition 2024.
- Check formatting with `cargo fmt --check` and lint all targets with `cargo clippy --all-targets`.
- Run the full suite with `cargo test`.
- Run one integration test with `cargo test --test dedup <test_name>`; run a CLI unit test with `cargo test --bin methdup <test_name>`.
- Tests build BAM fixtures with `noodles` in temporary directories; they do not require `samtools`, fixture downloads, or external services.

## Architecture

- `src/main.rs` owns argument parsing, stderr reporting, and exit codes; the public API is `deduplicate` re-exported from `src/lib.rs` and implemented in `src/processor.rs`.
- `src/alignment.rs` owns pair eligibility and exact duplicate-key construction; duplicate keys include both ends' reference, start, strand, and CIGAR plus the pair's library resolved from its `RG` tag via the header's `@RG.LB` (falling back to a shared `Unknown Library`). `src/pairing.rs` assembles reciprocal read1/read2 records by query name and library; `src/processor.rs` groups pairs and preserves coordinate order while writing.
- `Pairing::push` must emit each `Read` before its corresponding `Pair` or `Release`; the processor depends on the record already being in its ordered buffer.
- Duplicate groups retain the greatest combined base-quality pair, break ties by earliest input order, and recompute duplicate flags rather than preserving existing `0x400` flags.

## Safety Contracts

- Preserve the input-order buffer and coordinate-frontier finalization; writing resolved records directly can reorder the BAM or finalize a duplicate group too early.
- Output is written to a hidden sibling temporary BAM and renamed only after the writer finishes. Any validation, sorting, cache, or I/O error must leave the requested output unpublished.
- The header must declare coordinate sorting and records are checked independently for actual coordinate order. Ineligible, incomplete, and cross-reference records pass through unchanged.
- CLI processing failures print to stderr and exit 2; successful runs print the counter summary to stderr and exit 0.
