// Links module - derives edges between nodes in the variation graph
//
// This module determines the graph topology by finding which nodes connect
// to which other nodes based on the input sequences.

use std::io;
use std::sync::{Arc, RwLock};

use rayon::prelude::*;

use crate::intervaltree::AdaptiveTree;
use crate::intervaltree::IntervalTree;
use crate::pos::{incr_pos_by, is_rev, make_pos_t, offset, PosT};
use crate::seqindex::SeqIndex;

/// Represents a ranked/select-capable bitvector
/// This is a simplified interface for the SDSL sd_vector used in C++
#[derive(Clone)]
pub struct RankSelectBitVector {
    /// Positions where bits are set (sorted)
    pub set_positions: Vec<usize>,
    /// Total size of the bitvector
    pub size: usize,
}

impl RankSelectBitVector {
    /// Create from a raw bitvector
    pub fn from_bitvec(bits: &[bool]) -> Self {
        let set_positions: Vec<usize> = bits
            .iter()
            .enumerate()
            .filter_map(|(i, &b)| if b { Some(i) } else { None })
            .collect();

        RankSelectBitVector {
            set_positions,
            size: bits.len(),
        }
    }

    /// Rank: count number of 1s up to (but not including) position i
    pub fn rank(&self, i: usize) -> usize {
        self.set_positions.binary_search(&i).unwrap_or_else(|x| x)
    }

    /// Select: find position of the i-th 1 (1-indexed like SDSL)
    pub fn select(&self, i: usize) -> Option<usize> {
        if i == 0 || i > self.set_positions.len() {
            None
        } else {
            Some(self.set_positions[i - 1])
        }
    }

    /// Get the size of the bitvector
    pub fn size(&self) -> usize {
        self.size
    }

    /// Access: check if bit at position i is set
    pub fn access(&self, i: usize) -> bool {
        self.set_positions.binary_search(&i).is_ok()
    }
}

/// A simple set for storing links (pairs of positions)
/// In C++, this uses mmmulti::set, but we'll use a Vec for now
pub struct LinkSet {
    links: Vec<(PosT, PosT)>,
}

impl Default for LinkSet {
    fn default() -> Self {
        Self::new()
    }
}

impl LinkSet {
    pub fn new() -> Self {
        LinkSet { links: Vec::new() }
    }

    pub fn append(&mut self, link: (PosT, PosT)) {
        self.links.push(link);
    }

    pub fn links(&self) -> &[(PosT, PosT)] {
        &self.links
    }

    pub fn sort(&mut self) {
        self.links.par_sort_unstable();
    }

    /// Deduplicate links after sorting
    pub fn dedup(&mut self) {
        self.links.dedup();
    }

    pub fn len(&self) -> usize {
        self.links.len()
    }

    pub fn is_empty(&self) -> bool {
        self.links.is_empty()
    }
}

/// Derive links between nodes in the variation graph
///
/// For each node (marked in seq_id_cbv), determines which other nodes it
/// connects to by following the input sequences through the graph.
///
/// # Arguments
/// * `seqidx` - Sequence index
/// * `node_iitree` - Interval tree mapping graph positions to input positions
/// * `path_iitree` - Interval tree mapping input positions to graph positions
/// * `seq_id_cbv` - Ranked/select bitvector marking node boundaries
/// * `_num_threads` - Number of threads for parallel processing
///
/// # Returns
/// A LinkSet containing all edges in the graph as (from_node, to_node) pairs
pub fn derive_links(
    seqidx: Arc<SeqIndex>,
    node_iitree: Arc<RwLock<AdaptiveTree<u64, PosT>>>,
    path_iitree: Arc<RwLock<AdaptiveTree<u64, PosT>>>,
    seq_id_cbv: &RankSelectBitVector,
    _num_threads: usize,
) -> io::Result<LinkSet> {
    let n_nodes = seq_id_cbv.rank(seq_id_cbv.size() - 1);

    // Collect links in parallel
    // PERFORMANCE FIX: Use RwLock::read() instead of Mutex for concurrent read access
    let links: Vec<Vec<(PosT, PosT)>> = (1..=n_nodes)
        .into_par_iter()
        .map(|id| {
            let mut local_links = Vec::new();

            let node_start_in_s = match seq_id_cbv.select(id) {
                Some(pos) => pos as u64,
                None => return local_links,
            };
            let node_end_in_s = match seq_id_cbv.select(id + 1) {
                Some(pos) => pos as u64,
                None => return local_links,
            };

            // Find overlaps in node_iitree
            if let Ok(node_guard) = node_iitree.read() {
                node_guard
                    .overlap(node_start_in_s, node_end_in_s, |_idx, ovlp_start_in_s, ovlp_end_in_s, base_pos_start_in_q| {
                        let ovlp_length = ovlp_end_in_s - ovlp_start_in_s;
                        let mut pos_start_in_q = base_pos_start_in_q;
                        let mut pos_end_in_q = pos_start_in_q;

                        incr_pos_by(&mut pos_end_in_q, (ovlp_length - (ovlp_end_in_s - node_end_in_s)) as usize);
                        incr_pos_by(&mut pos_start_in_q, (node_start_in_s - ovlp_start_in_s) as usize);

                        let curr_step_is_rev = is_rev(pos_start_in_q);
                        let start_in_q = if curr_step_is_rev {
                            offset(pos_end_in_q) + 1
                        } else {
                            offset(pos_start_in_q)
                        };
                        let end_in_q = if curr_step_is_rev {
                            offset(pos_start_in_q) + 1
                        } else {
                            offset(pos_end_in_q)
                        };

                        let seq_id = match seqidx.seq_id_at(start_in_q) {
                            Some(id) => id as u64,
                            None => return,
                        };
                        let id_end = match seqidx.seq_id_at(end_in_q - 1) {
                            Some(id) => id as u64,
                            None => return,
                        };

                        if seq_id != id_end {
                            return;
                        }

                        let seq_start = match seqidx.nth_seq_offset(seq_id as usize) {
                            Some(offset) => offset,
                            None => return,
                        };
                        let seq_len = match seqidx.nth_seq_length(seq_id as usize) {
                            Some(len) => len,
                            None => return,
                        };
                        let seq_end = seq_start + seq_len;

                        // Only consider cases within sequence boundaries
                        if end_in_q < seq_end {
                            if let Ok(path_guard) = path_iitree.read() {
                                path_guard
                                    .overlap(end_in_q, end_in_q + 1, |_idx, ovlp_start_in_q, _ovlp_end_in_q, base_pos_start_in_s| {
                                        let mut pos_start_in_s = base_pos_start_in_s;
                                        if ovlp_start_in_q < end_in_q {
                                            incr_pos_by(&mut pos_start_in_s, (end_in_q - ovlp_start_in_q) as usize);
                                        }

                                        // Find which node we're in and record a link
                                        let next_id = seq_id_cbv.rank(offset(pos_start_in_s) as usize + 1);
                                        let next_step_is_rev = is_rev(pos_start_in_s);

                                        local_links.push((
                                            make_pos_t(id as u64, curr_step_is_rev),
                                            make_pos_t(next_id as u64, next_step_is_rev),
                                        ));
                                    })
                                    .ok();
                            }
                        }
                    })
                    .ok();
            }

            local_links
        })
        .collect();

    // Flatten and create LinkSet
    let mut link_set = LinkSet::new();
    for local_links in links {
        for link in local_links {
            link_set.append(link);
        }
    }

    // Sort and deduplicate the links
    link_set.sort();
    link_set.dedup();

    Ok(link_set)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_rank_select_bitvector() {
        let bits = vec![
            true, false, false, true, false, true, true, false, false, true,
        ];
        let bv = RankSelectBitVector::from_bitvec(&bits);

        // Test rank
        assert_eq!(bv.rank(0), 0); // No 1s before position 0
        assert_eq!(bv.rank(1), 1); // One 1 before position 1
        assert_eq!(bv.rank(4), 2); // Two 1s before position 4

        // Test select (1-indexed)
        assert_eq!(bv.select(1), Some(0)); // First 1 is at position 0
        assert_eq!(bv.select(2), Some(3)); // Second 1 is at position 3
        assert_eq!(bv.select(3), Some(5)); // Third 1 is at position 5
        assert_eq!(bv.select(0), None); // Invalid index
    }

    #[test]
    fn test_link_set() {
        let mut link_set = LinkSet::new();
        link_set.append((make_pos_t(1, false), make_pos_t(2, false)));
        link_set.append((make_pos_t(2, false), make_pos_t(3, true)));

        assert_eq!(link_set.len(), 2);

        link_set.sort();
        assert_eq!(link_set.len(), 2);
    }
}
