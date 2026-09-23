//! Command-line interface for methdup.

use std::{path::PathBuf, process::ExitCode};

use clap::Parser;
use methdup::{Counters, deduplicate};

#[derive(Debug, Parser)]
#[command(
    name = "methdup",
    version,
    about = "One-pass duplicate marking for coordinate-sorted BAM files"
)]
struct Cli {
    /// Input BAM file.
    input_bam: PathBuf,
    /// Output BAM file.
    output_bam: PathBuf,
    /// Omit losing pairs instead of writing them with SAM flag 0x400.
    #[arg(long)]
    remove_duplicates: bool,
    /// Maximum buffered records before failing safely.
    #[arg(long, default_value_t = 1_000_000)]
    max_cache_records: u64,
}

fn main() -> ExitCode {
    run(Cli::parse())
}

/// Runs deduplication and reports results, returning the process status.
fn run(cli: Cli) -> ExitCode {
    let counters = match deduplicate(
        &cli.input_bam,
        &cli.output_bam,
        cli.remove_duplicates,
        cli.max_cache_records,
    ) {
        Ok(counters) => counters,
        Err(error) => {
            eprintln!("methdup: error: {error}");
            return ExitCode::from(2);
        }
    };
    eprintln!("{}", format_summary(&counters));
    ExitCode::SUCCESS
}

/// Formats all processing counters as the single success-line CLI report.
pub fn format_summary(counters: &Counters) -> String {
    format!(
        "methdup: records={} eligible_pairs={} duplicate_pairs={} \
         flagged_records={} removed_records={} passthrough_records={} \
         peak_cache_records={}",
        counters.total_records,
        counters.eligible_pairs,
        counters.duplicate_pairs,
        counters.flagged_records,
        counters.removed_records,
        counters.passthrough_records,
        counters.peak_cache_records,
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_defaults() {
        let cli = Cli::try_parse_from(["methdup", "in.bam", "out.bam"]).expect("parses");
        assert!(!cli.remove_duplicates);
        assert_eq!(cli.max_cache_records, 1_000_000);
    }

    #[test]
    fn invalid_cache_limit_is_reported_as_failure() {
        let cli = Cli::try_parse_from(["methdup", "in.bam", "out.bam", "--max-cache-records", "0"])
            .expect("parses");
        // Exit status 2 mirrors the Python CLI contract.
        assert_eq!(run(cli), ExitCode::from(2));
    }
}
