# methdup

`methdup` deduplicates complete, same-reference proper pairs from a
coordinate-sorted BAM in one pass. By default it retains every record and sets
the SAM duplicate bit (`0x400`) on the lower-quality pair in each duplicate
group.

## Install

```sh
cargo install --path .
```

## Usage

```sh
methdup input.bam output.bam
methdup input.bam deduplicated.bam --remove-duplicates
```

On success it prints a summary to stderr and exits 0:

```text
methdup: records=4 eligible_pairs=2 duplicate_pairs=1 flagged_records=2 removed_records=0 passthrough_records=0 peak_cache_records=4
```

Duplicate identity includes both read ends' start coordinates, orientations,
and CIGAR strings, plus the pair's library resolved from its read group's `LB`
field (`RG` tag → `@RG`). Records without a resolvable library all share an
`Unknown Library` group, so they are compared against each other but never
against records from a named library. The representative is the pair with the
largest combined base-quality sum; input order breaks ties.

The command preserves coordinate order using a sliding cache. Set
`--max-cache-records N` to bound that cache (default: 1,000,000); it fails
safely if the bound is exceeded. It requires a BAM whose header declares
coordinate sorting and verifies that the records actually obey that ordering.

Unmapped, secondary, supplementary, incomplete, and inter-reference records
are passed through unchanged. No BAM index is created.
