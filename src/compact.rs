// Compact module - compresses non-bifurcating regions into nodes
//
// This module marks node boundaries in the variation graph by identifying
// positions where the graph structure changes (bifurcations, joins, etc.)

use std::io;
use std::sync::{Arc, RwLock};

use bitvec::prelude::*;
use rayon::prelude::*;

use std::sync::atomic::{AtomicU64, Ordering};

use crate::intervaltree::AdaptiveTree;
use crate::intervaltree::IntervalTree;
use crate::pos::{incr_pos_by, is_rev, offset, PosT};
use crate::seqindex::SeqIndex;

/// Atomic bitvector for thread-safe bit marking
struct AtomicBitVec {
    words: Vec<AtomicU64>,
    len: usize,
}

impl AtomicBitVec {
    fn new(size: usize) -> Self {
        let n_words = (size + 63) / 64;
        let words = (0..n_words).map(|_| AtomicU64::new(0)).collect();
        AtomicBitVec { words, len: size }
    }

    /// Atomically set a bit using fetch_or
    fn set(&self, index: usize) {
        let word_index = index / 64;
        let bit_index = index % 64;
        let mask = 1u64 << bit_index;
        self.words[word_index].fetch_or(mask, Ordering::Relaxed);
    }

    /// Get value of a bit
    #[allow(dead_code)]
    fn get(&self, index: usize) -> bool {
        let word_index = index / 64;
        let bit_index = index % 64;
        let mask = 1u64 << bit_index;
        (self.words[word_index].load(Ordering::Relaxed) & mask) != 0
    }

    /// Iterate over set bits
    fn iter_ones(&self) -> impl Iterator<Item = usize> + '_ {
        let len = self.len;
        self.words
            .iter()
            .enumerate()
            .flat_map(move |(word_idx, word)| {
                let w = word.load(Ordering::Relaxed);
                let base = word_idx * 64;
                (0..64).filter_map(move |bit| {
                    let idx = base + bit;
                    if idx < len && (w >> bit) & 1 == 1 {
                        Some(idx)
                    } else {
                        None
                    }
                })
            })
    }
}

/// Compact nodes by marking node boundaries in the graph
///
/// For each position in the input sequences, determines where it maps in the
/// graph sequence and marks the start and end of each mapping as a node boundary.
/// This compresses non-bifurcating linear regions into single nodes.
///
/// # Arguments
/// * `seqidx` - Sequence index
/// * `graph_size` - Size of the graph sequence
/// * `node_iitree` - Interval tree mapping graph positions to input positions
/// * `path_iitree` - Interval tree mapping input positions to graph positions
/// * `seq_id_bv` - Output bitvector marking node boundaries (will be modified)
/// * `num_threads` - Number of threads for parallel processing
pub fn compact_nodes(
    seqidx: Arc<SeqIndex>,
    graph_size: usize,
    _node_iitree: Arc<RwLock<AdaptiveTree<u64, PosT>>>,
    path_iitree: Arc<RwLock<AdaptiveTree<u64, PosT>>>,
    seq_id_bv: &mut BitVec<u64, Lsb0>,
    _num_threads: usize,
) -> io::Result<()> {
    // Create atomic bitvector for thread-safe marking
    let seq_id_abv = Arc::new(AtomicBitVec::new(seq_id_bv.len()));

    // Mark first position
    seq_id_abv.set(0);

    let num_seqs = seqidx.n_seqs();

    // Process each sequence in parallel
    (1..=num_seqs).into_par_iter().for_each(|i| {
        let j_start = match seqidx.nth_seq_offset(i) {
            Some(offset) => offset,
            None => return,
        };
        let seq_len = match seqidx.nth_seq_length(i) {
            Some(len) => len,
            None => return,
        };
        let k = j_start + seq_len;

        let mut j = j_start;
        while j < k {
            let mut overlap_count = 0u64;
            let mut ovlp_start_in_q = 0u64;
            let mut ovlp_end_in_q = 0u64;
            let mut pos_start_in_s = 0u64;

            // Find overlaps for this position
            if let Ok(path_guard) = path_iitree.read() {
                path_guard
                    .overlap(j, j + 1, |_idx, start, end, pos| {
                        overlap_count += 1;
                        ovlp_start_in_q = start;
                        ovlp_end_in_q = end;
                        pos_start_in_s = pos;
                    })
                    .ok();
            }

            // Each input base should map to exactly one place in the graph
            if overlap_count != 1 {
                let seq_name = seqidx
                    .nth_name(i)
                    .unwrap_or_else(|| "<unknown>".to_string());
                eprintln!(
                    "[compact] error: found {overlap_count} overlaps for seq {seq_name} idx {i} at j={j} of {k}"
                );
                if let Ok(path_guard) = path_iitree.read() {
                    path_guard
                        .overlap(j, j + 1, |_idx, start, end, pos| {
                            eprintln!(
                                "ovlp_start_in_q = {} ovlp_end_in_q = {} pos_start_in_s = {}+",
                                start,
                                end,
                                offset(pos)
                            );
                        })
                        .ok();
                }
                panic!("Overlap count mismatch");
            }

            let match_is_rev = is_rev(pos_start_in_s);
            let mut pos_end_in_s = pos_start_in_s;

            if !match_is_rev {
                // Forward match
                incr_pos_by(
                    &mut pos_end_in_s,
                    (ovlp_end_in_q - ovlp_start_in_q) as usize,
                );
                seq_id_abv.set(offset(pos_start_in_s) as usize);
                seq_id_abv.set(offset(pos_end_in_s) as usize);
            } else {
                // Reverse match
                incr_pos_by(
                    &mut pos_end_in_s,
                    (ovlp_end_in_q - ovlp_start_in_q - 1) as usize,
                );
                seq_id_abv.set(offset(pos_end_in_s) as usize);
                seq_id_abv.set((offset(pos_start_in_s) + 1) as usize);
            }

            j = ovlp_end_in_q;
        }
    });

    // Mark end of graph
    seq_id_abv.set(graph_size);

    // Copy atomic bitvector to output
    for pos in seq_id_abv.iter_ones() {
        seq_id_bv.set(pos, true);
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_atomic_bitvec() {
        let bv = AtomicBitVec::new(100);
        bv.set(5);
        bv.set(10);
        bv.set(50);

        assert!(bv.get(5));
        assert!(bv.get(10));
        assert!(bv.get(50));
        assert!(!bv.get(0));
        assert!(!bv.get(99));
    }

    #[test]
    fn test_atomic_bitvec_concurrent_same_word() {
        // Stress test: many threads setting different bits in the same u64 word.
        // The old read_volatile/write_volatile implementation would lose bits here.
        use std::sync::Arc;
        use std::thread;

        for _ in 0..100 {
            let bv = Arc::new(AtomicBitVec::new(64));
            let mut handles = Vec::new();
            for bit in 0..64 {
                let bv = Arc::clone(&bv);
                handles.push(thread::spawn(move || {
                    bv.set(bit);
                }));
            }
            for h in handles {
                h.join().unwrap();
            }
            for bit in 0..64 {
                assert!(bv.get(bit), "bit {} was lost due to race condition", bit);
            }
        }
    }
}
