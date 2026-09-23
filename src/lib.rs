//! Public API for the methdup BAM deduplicator.

pub mod alignment;
pub mod models;
pub mod pairing;
pub mod processor;

pub use models::{
    BufferedRecord, BufferedRecordInner, Counters, DeduplicationError, Pair, PairKey,
};
pub use processor::deduplicate;
