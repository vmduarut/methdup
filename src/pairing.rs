//! Streaming assembly of primary paired-end alignments.

use std::{cell::RefCell, collections::HashMap, rc::Rc};

use bstr::BString;
use noodles::sam::{Header, alignment::RecordBuf};

use crate::{
    alignment::{
        build_pair_key, have_reciprocal_mate_coordinates, is_eligible_pair_end, record_library,
    },
    models::{BufferedRecord, BufferedRecordInner, PairKey},
};

/// An input alignment that must be added to the ordered output buffer.
#[derive(Debug)]
pub enum PairingEvent {
    /// A read to append to the ordered output buffer.
    Read {
        item: BufferedRecord,
        passthrough: bool,
    },
    /// A complete reciprocal pair whose two reads were already emitted.
    Pair {
        first: BufferedRecord,
        second: BufferedRecord,
        key: PairKey,
        rightmost_start: i32,
        input_order: u64,
    },
    /// An eligible record which cannot form a valid pair and may be written.
    Release { item: BufferedRecord },
}

/// Streams records into read, completed-pair, and safe-release events.
///
/// Every input record first produces a read event for ordered buffering.
/// Eligible mates are retained by query name and library until a reciprocal
/// read1/read2 pair completes, or released when coordinate order proves a mate
/// cannot arrive. Invalid, incomplete, and non-candidate records are released
/// for unchanged output. Keying by both query name and library keeps
/// same-name records from different libraries from ever crossing pairs.
#[derive(Debug)]
pub struct Pairing<'a> {
    header: &'a Header,
    pending: HashMap<(BString, BString), [Option<BufferedRecord>; 2]>,
    invalid_names: std::collections::HashSet<(BString, BString)>,
    record_sequence: u64,
    current_reference: Option<i32>,
    current_start: Option<i32>,
}

impl<'a> Pairing<'a> {
    /// Creates an empty pairing stream bound to a BAM header.
    pub fn new(header: &'a Header) -> Self {
        Self {
            header,
            pending: HashMap::new(),
            invalid_names: std::collections::HashSet::new(),
            record_sequence: 0,
            current_reference: None,
            current_start: None,
        }
    }

    /// Feeds one input record, returning the events it produces in order.
    pub fn push(&mut self, record: RecordBuf) -> Vec<PairingEvent> {
        let mut events = Vec::new();
        let unmapped = record.flags().is_unmapped();
        let candidate = is_eligible_pair_end(&record);

        let reference = if unmapped {
            None
        } else {
            Some(crate::alignment::reference_id(&record))
        };
        let start = if unmapped {
            None
        } else {
            Some(crate::alignment::reference_start(&record))
        };

        if unmapped {
            self.append_expired(&mut events, (None, None));
            self.invalid_names.clear();
            self.current_reference = None;
            self.current_start = None;
        } else if self.current_reference != reference {
            self.append_expired(&mut events, (None, None));
            self.invalid_names.clear();
            self.current_reference = reference;
            self.current_start = start;
        } else if self.current_start != start {
            self.append_expired(&mut events, (reference, start));
            self.current_start = start;
        }

        self.record_sequence += 1;
        let item: BufferedRecord = Rc::new(RefCell::new(BufferedRecordInner {
            record,
            blocked: candidate,
            ordinal: self.record_sequence,
            drop: false,
        }));

        events.push(PairingEvent::Read {
            item: Rc::clone(&item),
            passthrough: !candidate,
        });

        if !candidate {
            return events;
        }

        let name: BString = item.borrow().record.name().expect("eligible name").into();
        let library = record_library(self.header, &item.borrow().record);
        let key = (name, library);
        if self.invalid_names.contains(&key) {
            events.push(PairingEvent::Release { item });
            return events;
        }

        let end_index = if item.borrow().record.flags().is_first_segment() {
            0
        } else {
            1
        };
        let other_index = 1 - end_index;

        if !self.pending.contains_key(&key) {
            self.pending.insert(key.clone(), [None, None]);
        }
        let ends = self.pending.get_mut(&key).expect("slot exists");

        if ends[end_index].is_some() {
            let previous = ends[end_index].take().expect("occupied");
            if ends[other_index].is_none() {
                self.pending.remove(&key);
            }
            self.invalid_names.insert(key);
            events.push(PairingEvent::Release { item: previous });
            events.push(PairingEvent::Release { item });
            return events;
        }

        if let Some(other) = ends[other_index].take() {
            self.pending.remove(&key);
            let (read1_item, read2_item) = if end_index == 0 {
                (&item, &other)
            } else {
                (&other, &item)
            };
            if !have_reciprocal_mate_coordinates(
                &read1_item.borrow().record,
                &read2_item.borrow().record,
            ) {
                events.push(PairingEvent::Release {
                    item: Rc::clone(read1_item),
                });
                events.push(PairingEvent::Release {
                    item: Rc::clone(read2_item),
                });
                return events;
            }
            let read1_rec = &read1_item.borrow().record;
            let read2_rec = &read2_item.borrow().record;
            events.push(PairingEvent::Pair {
                first: Rc::clone(read1_item),
                second: Rc::clone(read2_item),
                key: build_pair_key(read1_rec, read2_rec, self.header),
                rightmost_start: read1_rec
                    .alignment_start()
                    .zip(read2_rec.alignment_start())
                    .map_or(-1, |(a, b)| {
                        i32::try_from(a.get().max(b.get()) - 1).unwrap_or(i32::MAX)
                    }),
                input_order: read1_item.borrow().ordinal.min(read2_item.borrow().ordinal),
            });
            return events;
        }

        let next_reference_start = crate::alignment::mate_start(&item.borrow().record);
        if let Some(current_start) = self.current_start
            && next_reference_start < current_start
        {
            self.pending.remove(&key);
            events.push(PairingEvent::Release { item });
            return events;
        }

        ends[end_index] = Some(item);
        events
    }

    /// Releases every remaining pending record; call once at end of input.
    pub fn finish(&mut self) -> Vec<PairingEvent> {
        let mut events = Vec::new();
        self.append_expired(&mut events, (None, None));
        events
    }

    fn append_expired(
        &mut self,
        events: &mut Vec<PairingEvent>,
        (next_reference, next_start): (Option<i32>, Option<i32>),
    ) {
        let mut released_keys = Vec::new();
        for (key, ends) in self.pending.iter_mut() {
            for end in ends.iter_mut() {
                if let Some(item) = end {
                    let passed_mate = {
                        let record = &item.borrow().record;
                        next_reference.is_none()
                            || crate::alignment::reference_id(record)
                                != next_reference.unwrap_or(-1)
                            || (next_start.is_some()
                                && crate::alignment::mate_start(record) < next_start.unwrap_or(-1))
                    };
                    if passed_mate {
                        events.push(PairingEvent::Release {
                            item: Rc::clone(item),
                        });
                        *end = None;
                    }
                }
            }
            if ends.iter().all(Option::is_none) {
                released_keys.push(key.clone());
            }
        }
        for key in released_keys {
            self.pending.remove(&key);
        }
    }
}
