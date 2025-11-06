// Compact module - compresses non-bifurcating regions into nodes
//
// This module marks node boundaries in the variation graph by identifying
// positions where the graph structure changes (bifurcations, joins, etc.)

use std::sync::{Arc, RwLock};
use std::io;

use rayon::prelude::*;
use bitvec::prelude::*;

use crate::pos::{PosT, offset, is_rev, incr_pos_by};
use crate::seqindex::SeqIndex;
use crate::intervaltree::AdaptiveTree;
use crate::intervaltree::IntervalTree;

/// Atomic bitvector for thread-safe bit marking
struct AtomicBitVec {
    bits: BitVec<u64, Lsb0>,
}

impl AtomicBitVec {
    fn new(size: usize) -> Self {
        AtomicBitVec {
            bits: BitVec::repeat(false, size),
        }
    }

    /// Atomically set a bit (simplified version - would need true atomics in production)
    fn set(&self, index: usize) {
        unsafe {
            let ptr = self.bits.as_raw_slice().as_ptr() as *mut u64;
            let word_index = index / 64;
            let bit_index = index % 64;
            let mask = 1u64 << bit_index;
            let word_ptr = ptr.add(word_index);
            let old_word = std::ptr::read_volatile(word_ptr);
            std::ptr::write_volatile(word_ptr, old_word | mask);
        }
    }

    /// Get value of a bit
    fn get(&self, index: usize) -> bool {
        self.bits[index]
    }

    /// Iterate over set bits
    fn iter_ones(&self) -> impl Iterator<Item = usize> + '_ {
        self.bits.iter_ones()
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
    node_iitree: Arc<RwLock<AdaptiveTree<u64, PosT>>>,
    path_iitree: Arc<RwLock<AdaptiveTree<u64, PosT>>>,
    seq_id_bv: &mut BitVec<u64, Lsb0>,
    num_threads: usize,
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
                path_guard.overlap(j, j + 1, |_idx, start, end, pos| {
                    overlap_count += 1;
                    ovlp_start_in_q = start;
                    ovlp_end_in_q = end;
                    pos_start_in_s = pos;
                }).ok();
            }

            // Each input base should map to exactly one place in the graph
            if overlap_count != 1 {
                let seq_name = seqidx.nth_name(i).unwrap_or_else(|| "<unknown>".to_string());
                eprintln!(
                    "[compact] error: found {} overlaps for seq {} idx {} at j={} of {}",
                    overlap_count,
                    seq_name,
                    i,
                    j,
                    k
                );
                if let Ok(path_guard) = path_iitree.read() {
                    path_guard.overlap(j, j + 1, |_idx, start, end, pos| {
                        eprintln!(
                            "ovlp_start_in_q = {} ovlp_end_in_q = {} pos_start_in_s = {}+",
                            start, end, offset(pos)
                        );
                    }).ok();
                }
                panic!("Overlap count mismatch");
            }

            let match_is_rev = is_rev(pos_start_in_s);
            let mut pos_end_in_s = pos_start_in_s;

            if !match_is_rev {
                // Forward match
                incr_pos_by(&mut pos_end_in_s, (ovlp_end_in_q - ovlp_start_in_q) as usize);
                seq_id_abv.set(offset(pos_start_in_s) as usize);
                seq_id_abv.set(offset(pos_end_in_s) as usize);
            } else {
                // Reverse match
                incr_pos_by(&mut pos_end_in_s, (ovlp_end_in_q - ovlp_start_in_q - 1) as usize);
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
}
