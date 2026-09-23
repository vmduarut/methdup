//! Streaming, coordinate-order-preserving BAM deduplication.

use std::{
    cmp::Reverse,
    collections::{BinaryHeap, HashMap, VecDeque},
    fs::File,
    path::Path,
    rc::Rc,
};

use noodles::{
    bam::{self},
    bgzf::{self, io::Writer as BgzfWriter},
    sam::{
        self,
        alignment::{RecordBuf, io::Write as _},
        header::record::value::map::header::tag::SORT_ORDER,
    },
};

use crate::{
    alignment::{add_program_record, set_duplicate},
    models::{BufferedRecord, Counters, DeduplicationError, Pair, PairKey},
    pairing::{Pairing, PairingEvent},
};

/// Validates BAM suffixes and rejects an output equal to the input path.
pub fn validate_paths(input: &Path, output: &Path) -> Result<(), DeduplicationError> {
    let has_bam_extension = |path: &Path| {
        path.extension()
            .and_then(|extension| extension.to_str())
            .is_some_and(|extension| extension.eq_ignore_ascii_case("bam"))
    };
    if !has_bam_extension(input) || !has_bam_extension(output) {
        return Err(DeduplicationError::InvalidPaths);
    }
    let canonical = |path: &Path| -> std::path::PathBuf {
        path.canonicalize().unwrap_or_else(|_| path.to_path_buf())
    };
    if canonical(input) == canonical(output) {
        return Err(DeduplicationError::SamePath);
    }
    Ok(())
}

/// Opens a BAM reader only when its header declares coordinate sorting.
pub fn open_coordinate_sorted_bam(
    path: &Path,
) -> Result<(bam::io::Reader<bgzf::io::Reader<File>>, sam::Header), DeduplicationError> {
    let file = File::open(path)?;
    let mut reader = bam::io::Reader::new(file);
    let header = reader.read_header()?;
    let sort_order = header
        .header()
        .and_then(|header| header.other_fields().get(&SORT_ORDER));
    if sort_order.is_none_or(|value| value != "coordinate") {
        return Err(DeduplicationError::UncoordinateSorted);
    }
    Ok((reader, header))
}

/// A BAM output stream written to a temporary file and atomically published.
pub struct AtomicBamWriter {
    output: std::path::PathBuf,
    temp_path: Option<std::path::PathBuf>,
    writer: Option<bam::io::Writer<BgzfWriter<File>>>,
}

impl AtomicBamWriter {
    pub fn begin(output: &Path) -> Result<Self, DeduplicationError> {
        if let Some(parent) = output.parent()
            && !parent.as_os_str().is_empty()
        {
            std::fs::create_dir_all(parent)?;
        }
        let file_name = output
            .file_name()
            .and_then(|name| name.to_str())
            .unwrap_or("output");
        let mut attempt = 0u32;
        let temp_path = loop {
            let candidate = output.with_file_name(format!(
                ".{file_name}.{}.{}.tmp.bam",
                std::process::id(),
                attempt
            ));
            match File::options()
                .write(true)
                .create_new(true)
                .open(&candidate)
            {
                Ok(_) => break candidate,
                Err(error) if error.kind() == std::io::ErrorKind::AlreadyExists => attempt += 1,
                Err(error) => return Err(error.into()),
            }
        };
        let file = File::options().write(true).open(&temp_path)?;
        Ok(Self {
            output: output.to_path_buf(),
            temp_path: Some(temp_path),
            writer: Some(bam::io::Writer::new(file)),
        })
    }

    pub fn write_header(&mut self, header: &sam::Header) -> Result<(), DeduplicationError> {
        self.writer
            .as_mut()
            .expect("uncommitted writer")
            .write_header(header)?;
        Ok(())
    }

    pub fn write_record(
        &mut self,
        header: &sam::Header,
        record: &RecordBuf,
    ) -> Result<(), DeduplicationError> {
        self.writer
            .as_mut()
            .expect("uncommitted writer")
            .write_alignment_record(header, record)?;
        Ok(())
    }

    pub fn commit(mut self) -> Result<(), DeduplicationError> {
        let mut writer = self.writer.take().expect("uncommitted writer");
        writer.try_finish()?;
        drop(writer);
        let temp_path = self.temp_path.take().expect("uncommitted writer");
        std::fs::rename(&temp_path, &self.output)?;
        Ok(())
    }
}

impl Drop for AtomicBamWriter {
    fn drop(&mut self) {
        if let Some(path) = self.temp_path.take() {
            let _ = std::fs::remove_file(path);
        }
    }
}

/// Deduplicates a coordinate-sorted BAM and atomically publishes the result.
pub fn deduplicate(
    input: &Path,
    output: &Path,
    remove_duplicates: bool,
    max_cache_records: u64,
) -> Result<Counters, DeduplicationError> {
    if max_cache_records < 1 {
        return Err(DeduplicationError::InvalidCacheLimit);
    }
    validate_paths(input, output)?;

    let (mut reader, mut header) = open_coordinate_sorted_bam(input)?;
    add_program_record(&mut header)?;
    let mut atomic = AtomicBamWriter::begin(output)?;
    atomic.write_header(&header)?;

    let counters = {
        let mut processor =
            Processor::new(&mut atomic, &header, remove_duplicates, max_cache_records);
        processor.process(&mut reader)?;
        processor.counters
    };
    atomic.commit()?;
    Ok(counters)
}

/// Owns the mutable state required for one-pass ordered deduplication.
pub struct Processor<'a> {
    destination: &'a mut AtomicBamWriter,
    header: &'a sam::Header,
    counters: Counters,
    remove_duplicates: bool,
    max_cache_records: u64,
    buffer: VecDeque<BufferedRecord>,
    pairing: Pairing<'a>,
    groups: HashMap<PairKey, Vec<Pair>>,
    group_heap: BinaryHeap<Reverse<(i32, u64, PairKey)>>,
    group_sequence: u64,
    current_reference: Option<i32>,
    current_start: Option<i32>,
    last_coordinate: Option<(i32, i32)>,
    seen_unmapped: bool,
}

impl<'a> Processor<'a> {
    pub fn new(
        destination: &'a mut AtomicBamWriter,
        header: &'a sam::Header,
        remove_duplicates: bool,
        max_cache_records: u64,
    ) -> Self {
        Self {
            destination,
            header,
            counters: Counters::default(),
            remove_duplicates,
            max_cache_records,
            buffer: VecDeque::new(),
            pairing: Pairing::new(header),
            groups: HashMap::new(),
            group_heap: BinaryHeap::new(),
            group_sequence: 0,
            current_reference: None,
            current_start: None,
            last_coordinate: None,
            seen_unmapped: false,
        }
    }

    pub fn process(
        &mut self,
        reader: &mut bam::io::Reader<bgzf::io::Reader<File>>,
    ) -> Result<(), DeduplicationError> {
        for result in reader.record_bufs(self.header) {
            let record = result.map_err(DeduplicationError::Io)?;
            let events = self.pairing.push(record);
            self.handle_events(events)?;
        }
        let events = self.pairing.finish();
        self.handle_events(events)?;
        self.finalize_groups(None)?;
        self.flush_unblocked()?;
        if !self.buffer.is_empty() {
            return Err(DeduplicationError::BlockedAtEnd);
        }
        Ok(())
    }

    fn handle_events(&mut self, events: Vec<PairingEvent>) -> Result<(), DeduplicationError> {
        for event in events {
            match event {
                PairingEvent::Read { item, passthrough } => self.handle_read(item, passthrough)?,
                PairingEvent::Pair {
                    first,
                    second,
                    key,
                    rightmost_start,
                    input_order,
                } => self.handle_pair(first, second, key, rightmost_start, input_order),
                PairingEvent::Release { item } => self.handle_release(item)?,
            }
        }
        Ok(())
    }

    fn handle_read(
        &mut self,
        item: BufferedRecord,
        passthrough: bool,
    ) -> Result<(), DeduplicationError> {
        let record = item.borrow().record.clone();
        self.counters.total_records += 1;
        self.validate_and_advance_coordinate(&record)?;

        self.buffer.push_back(item);
        self.counters.peak_cache_records = self
            .counters
            .peak_cache_records
            .max(self.buffer.len() as u64);
        if self.buffer.len() as u64 > self.max_cache_records {
            return Err(DeduplicationError::CacheExceeded(self.max_cache_records));
        }
        if passthrough {
            self.counters.passthrough_records += 1;
        }
        self.flush_unblocked()?;
        Ok(())
    }

    fn handle_release(&mut self, item: BufferedRecord) -> Result<(), DeduplicationError> {
        item.borrow_mut().blocked = false;
        self.counters.passthrough_records += 1;
        self.flush_unblocked()?;
        Ok(())
    }

    fn handle_pair(
        &mut self,
        first: BufferedRecord,
        second: BufferedRecord,
        key: PairKey,
        rightmost_start: i32,
        input_order: u64,
    ) {
        self.counters.eligible_pairs += 1;
        if !self.groups.contains_key(&key) {
            self.group_sequence += 1;
            let sequence = self.group_sequence;
            self.group_heap
                .push(Reverse((rightmost_start, sequence, key.clone())));
        }
        self.groups.entry(key).or_default().push(Pair {
            first,
            second,
            rightmost_start,
            input_order,
        });
    }

    fn validate_and_advance_coordinate(
        &mut self,
        record: &RecordBuf,
    ) -> Result<(), DeduplicationError> {
        use crate::alignment::{reference_id, reference_start};
        if record.flags().is_unmapped() {
            self.seen_unmapped = true;
            self.advance_coordinate_frontier(None, None)?;
            return Ok(());
        }
        let coordinate = (reference_id(record), reference_start(record));
        if self.seen_unmapped || self.last_coordinate.is_some_and(|last| coordinate < last) {
            return Err(DeduplicationError::UnsortedRecords);
        }
        self.last_coordinate = Some(coordinate);
        self.advance_coordinate_frontier(Some(coordinate.0), Some(coordinate.1))?;
        Ok(())
    }

    fn advance_coordinate_frontier(
        &mut self,
        reference_id: Option<i32>,
        start: Option<i32>,
    ) -> Result<(), DeduplicationError> {
        if self.current_reference.is_none() {
            self.current_reference = reference_id;
            self.current_start = start;
            return Ok(());
        }
        if reference_id == self.current_reference && start == self.current_start {
            return Ok(());
        }
        if reference_id == self.current_reference && start.is_some() {
            self.finalize_groups(start)?;
        } else {
            self.finalize_groups(None)?;
        }
        self.flush_unblocked()?;
        self.current_reference = reference_id;
        self.current_start = start;
        Ok(())
    }

    fn finalize_groups(&mut self, before_start: Option<i32>) -> Result<(), DeduplicationError> {
        while self
            .group_heap
            .peek()
            .is_some_and(|Reverse((rightmost, _, _))| {
                before_start.is_none_or(|before| *rightmost < before)
            })
        {
            let Reverse((_, _, key)) = self.group_heap.pop().expect("peeked non-empty");
            if let Some(group) = self.groups.remove(&key) {
                self.finalize_duplicate_group(group);
            }
        }
        Ok(())
    }

    fn finalize_duplicate_group(&mut self, group: Vec<Pair>) {
        let winner = group
            .iter()
            .max_by_key(|pair| (pair.base_quality_sum(), Reverse(pair.input_order)))
            .expect("non-empty group");
        for pair in &group {
            let loser = !Rc::ptr_eq(&pair.first, &winner.first);
            set_duplicate(&mut pair.first.borrow_mut().record, loser);
            set_duplicate(&mut pair.second.borrow_mut().record, loser);
            if loser {
                self.counters.duplicate_pairs += 1;
                self.counters.flagged_records += 2;
                if self.remove_duplicates {
                    pair.first.borrow_mut().drop = true;
                    pair.second.borrow_mut().drop = true;
                }
            }
            pair.first.borrow_mut().blocked = false;
            pair.second.borrow_mut().blocked = false;
        }
    }

    fn flush_unblocked(&mut self) -> Result<(), DeduplicationError> {
        while let Some(front) = self.buffer.front() {
            if front.borrow().blocked {
                break;
            }
            let item = self.buffer.pop_front().expect("front exists");
            if item.borrow().drop {
                self.counters.removed_records += 1;
            } else {
                let record = item.borrow().record.clone();
                self.destination.write_record(self.header, &record)?;
            }
        }
        Ok(())
    }
}
