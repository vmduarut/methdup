//! Data structures shared by the duplicate-detection pipeline.

use std::cell::RefCell;
use std::rc::Rc;

use bstr::BString;
use noodles::sam::alignment::RecordBuf;

/// An immutable alignment identity used to group duplicate fragments.
///
/// The key keeps both ends in read1/read2 order and includes each end's
/// reference, start, strand, and CIGAR so alignments with different internal
/// structures are not collapsed into the same duplicate group.
#[derive(Clone, Debug, Eq, Hash, Ord, PartialEq, PartialOrd)]
pub struct PairKey {
    pub read1_reference: i32,
    pub read1_start: i32,
    pub read1_reverse: bool,
    pub read1_cigar: Option<String>,
    pub read2_reference: i32,
    pub read2_start: i32,
    pub read2_reverse: bool,
    pub read2_cigar: Option<String>,
    /// The resolved library (`@RG.LB`) both reads belong to.
    pub library: BString,
}

/// Processing counters reported by the CLI.
#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct Counters {
    pub total_records: u64,
    pub eligible_pairs: u64,
    pub duplicate_pairs: u64,
    pub flagged_records: u64,
    pub removed_records: u64,
    pub passthrough_records: u64,
    pub peak_cache_records: u64,
}

/// A buffered input alignment with its mutable resolution state.
///
/// Wrapped in a shared cell so the pairing stage, the ordered output buffer,
/// and duplicate groups can all observe and mutate the same record.
#[derive(Debug)]
pub struct BufferedRecordInner {
    pub record: RecordBuf,
    pub blocked: bool,
    pub ordinal: u64,
    pub drop: bool,
}

/// The shared, mutable identity of one buffered record.
pub type BufferedRecord = Rc<RefCell<BufferedRecordInner>>;

/// A complete reciprocal pair whose reads were already buffered.
#[derive(Debug)]
pub struct Pair {
    pub first: BufferedRecord,
    pub second: BufferedRecord,
    pub rightmost_start: i32,
    pub input_order: u64,
}

impl Pair {
    /// Returns the combined base-quality sum used to choose a group winner.
    ///
    /// Missing query-quality arrays contribute zero, allowing the processor to
    /// compare otherwise valid pairs deterministically.
    pub fn base_quality_sum(&self) -> u64 {
        self.first
            .borrow()
            .record
            .quality_scores()
            .iter()
            .map(u64::from)
            .sum::<u64>()
            + self
                .second
                .borrow()
                .record
                .quality_scores()
                .iter()
                .map(u64::from)
                .sum::<u64>()
    }
}

/// An error raised when input cannot safely be deduplicated in one pass.
#[derive(Debug, thiserror::Error)]
pub enum DeduplicationError {
    #[error("input and output paths must both end in .bam")]
    InvalidPaths,
    #[error("input and output paths must be different")]
    SamePath,
    #[error("input BAM header must declare coordinate sort order")]
    UncoordinateSorted,
    #[error("input BAM records are not coordinate sorted")]
    UnsortedRecords,
    #[error(
        "cache exceeded --max-cache-records ({0}); increase the limit or use inputs with shorter mate spans"
    )]
    CacheExceeded(u64),
    #[error("--max-cache-records must be at least 1")]
    InvalidCacheLimit,
    #[error("internal error: records remained blocked at end of input")]
    BlockedAtEnd,
    #[error("{0}")]
    Io(#[from] std::io::Error),
    #[error("{0}")]
    Other(String),
}
