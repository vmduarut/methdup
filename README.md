# methdup

`methdup` deduplicates complete, same-reference proper pairs from a
coordinate-sorted BAM in one pass. By default it retains every record and sets
the SAM duplicate bit (`0x400`) on the lower-quality pair in each duplicate
group.

```sh
methdup input.bam output.bam
methdup input.bam deduplicated.bam --remove-duplicates
```

Duplicate identity includes both read ends' start coordinates, orientations,
and CIGAR strings. The representative is the pair with the largest combined
base-quality sum; input order breaks ties.

The command preserves coordinate order using a sliding cache. Set
`--max-cache-records N` to bound that cache (default: 1,000,000); it fails
safely if the bound is exceeded. It requires a BAM whose header declares
coordinate sorting and verifies that the records actually obey that ordering.

Unmapped, secondary, supplementary, incomplete, and inter-reference records
are passed through unchanged. No BAM index is created.
