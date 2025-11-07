// Transclosure computation for variation graph construction
//
// This module computes transitive closures of aligned positions,
// identifying equivalence classes that form nodes in the variation graph.

use std::collections::HashMap;
use std::io;
use std::sync::atomic::Ordering;
use std::sync::{Arc, RwLock};
use std::thread;

use crate::dset64_asm::DisjointSetsAsm;
use crossbeam_queue::ArrayQueue;
use rayon::prelude::*;
use std::sync::atomic::AtomicU64;
use sucds::bit_vectors::{BitVector as SucdsBitVector, Rank, Rank9Sel};

use crate::intervaltree::AdaptiveTree;
use crate::intervaltree::IntervalTree;
use crate::pos::{decr_pos, decr_pos_by, incr_pos, incr_pos_by, is_rev, make_pos_t, offset, PosT};
use crate::seqindex::SeqIndex;

/// Thomas Wang's 64-bit integer hash function
///
/// In many implementations, std::hash is identity for integers,
/// which leads to performance issues. This provides better distribution.
#[inline]
pub fn wang_hash_64(mut key: u64) -> u64 {
    key = (!key).wrapping_add(key << 21); // key = (key << 21) - key - 1
    key = key ^ (key >> 24);
    key = key.wrapping_add(key << 3).wrapping_add(key << 8); // key * 265
    key = key ^ (key >> 14);
    key = key.wrapping_add(key << 2).wrapping_add(key << 4); // key * 21
    key = key ^ (key >> 28);
    key = key.wrapping_add(key << 31);
    key
}

/// Range in the graph sequence
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Range {
    pub begin: u64,
    pub end: u64,
}

impl Range {
    pub fn new(begin: u64, end: u64) -> Self {
        Range { begin, end }
    }
}

/// Match representing an alignment between positions
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Match {
    pub start: u64,
    pub end: u64,
    pub data: PosT,
}

impl Match {
    pub fn new(start: u64, end: u64, data: PosT) -> Self {
        Match { start, end, data }
    }

    pub fn length(&self) -> u64 {
        self.end - self.start
    }
}

/// Type aliases for atomic queues and bitvector
type RangeAtomicQueue = ArrayQueue<(PosT, u64)>;
type OverlapAtomicQueue = ArrayQueue<(Match, bool)>;

/// Atomic bitvector wrapper for thread-safe bit operations
/// Uses Vec<AtomicU64> with fetch_or for truly atomic bit operations
/// (compiles to lock-free CPU instructions like LOCK OR on x86)
#[derive(Debug)]
struct AtomicBitVec {
    data: Vec<AtomicU64>,
}

impl AtomicBitVec {
    fn new(size: usize) -> Self {
        let num_words = (size + 63) / 64;
        let mut data = Vec::with_capacity(num_words);
        for _ in 0..num_words {
            data.push(AtomicU64::new(0));
        }
        AtomicBitVec { data }
    }

    /// Atomically set a bit and return its previous value
    fn set(&self, index: usize, _value: bool, ordering: Ordering) -> bool {
        let word_index = index / 64;
        let bit_index = index % 64;
        let mask = 1u64 << bit_index;
        let prev = self.data[word_index].fetch_or(mask, ordering);
        (prev & mask) != 0
    }

    /// Get a bit value
    fn get(&self, index: usize, ordering: Ordering) -> bool {
        let word_index = index / 64;
        let bit_index = index % 64;
        let mask = 1u64 << bit_index;
        let word = self.data[word_index].load(ordering);
        (word & mask) != 0
    }
}

/// Extend a range in the buffer
///
/// Finds an existing range to extend or creates a new one.
/// Flushes ranges when hitting sequence boundaries.
pub fn extend_range(
    s_pos: u64,
    q_pos: PosT,
    range_buffer: &mut HashMap<PosT, Range>,
    seqidx: &SeqIndex,
    node_iitree: &mut AdaptiveTree<u64, PosT>,
    path_iitree: &mut AdaptiveTree<u64, PosT>,
) -> io::Result<()> {
    // Find position to add onto (must match position and orientation)
    let mut q_last_pos = q_pos;
    decr_pos(&mut q_last_pos);

    if let Some(found) = range_buffer.get(&q_last_pos).copied() {
        // Check if we're at a sequence boundary
        let at_boundary = if !is_rev(q_pos) {
            seqidx.seq_start(offset(q_pos))
        } else {
            seqidx.seq_start(offset(q_last_pos))
        };

        if at_boundary {
            // Flush the buffer we found (don't extend across node boundaries)
            flush_single_range(&found, q_last_pos, node_iitree, path_iitree)?;
            range_buffer.remove(&q_last_pos);
            range_buffer.insert(q_pos, Range::new(s_pos, s_pos + 1));
        } else if found.end == s_pos {
            // Extend the existing range
            range_buffer.remove(&q_last_pos);
            range_buffer.insert(q_pos, Range::new(found.begin, s_pos + 1));
        } else {
            // Store a new range
            range_buffer.insert(q_pos, Range::new(s_pos, s_pos + 1));
        }
    } else {
        // No existing range, create new one
        range_buffer.insert(q_pos, Range::new(s_pos, s_pos + 1));
    }

    Ok(())
}

/// Flush a single range to the iitrees
///
/// Given a range ending at match_end_pos_in_q, compute the match parameters
/// and add intervals to both node and path iitrees.
fn flush_single_range(
    range_in_s: &Range,
    match_end_pos_in_q: PosT,
    node_iitree: &mut AdaptiveTree<u64, PosT>,
    path_iitree: &mut AdaptiveTree<u64, PosT>,
) -> io::Result<()> {
    let is_rev_match = is_rev(match_end_pos_in_q);
    let match_length = range_in_s.end - range_in_s.begin;
    let match_start_in_s = range_in_s.begin;
    let match_end_in_s = range_in_s.end;

    let (match_pos_in_s, match_pos_in_q, match_start_in_q, match_end_in_q) = if !is_rev_match {
        // Forward match
        let match_end_in_q = offset(match_end_pos_in_q) + 1;
        let match_start_in_q = match_end_in_q - match_length;
        let match_pos_in_s = make_pos_t(match_start_in_s, false);
        let match_pos_in_q = make_pos_t(match_start_in_q, false);
        (
            match_pos_in_s,
            match_pos_in_q,
            match_start_in_q,
            match_end_in_q,
        )
    } else {
        // Reverse match
        let match_end_in_q = offset(match_end_pos_in_q);
        let mut match_end_pos_in_q_tmp = match_end_pos_in_q;
        decr_pos_by(&mut match_end_pos_in_q_tmp, match_length as usize);
        let match_pos_in_s = make_pos_t(match_end_in_s - 1, true);
        let match_pos_in_q = make_pos_t(offset(match_end_pos_in_q_tmp) - 1, true);
        let match_start_in_q = match_end_in_q;
        let match_end_in_q = offset(match_end_pos_in_q_tmp);
        (
            match_pos_in_s,
            match_pos_in_q,
            match_start_in_q,
            match_end_in_q,
        )
    };

    // Add to both iitrees
    node_iitree.add(match_start_in_s, match_end_in_s, match_pos_in_q);
    path_iitree.add(match_start_in_q, match_end_in_q, match_pos_in_s);

    Ok(())
}

/// Flush all ranges in the buffer that aren't at s_pos
pub fn flush_ranges(
    s_pos: u64,
    range_buffer: &mut HashMap<PosT, Range>,
    node_iitree: &mut AdaptiveTree<u64, PosT>,
    path_iitree: &mut AdaptiveTree<u64, PosT>,
) -> io::Result<()> {
    let to_flush: Vec<_> = range_buffer
        .iter()
        .filter(|(_, range)| range.end != s_pos)
        .map(|(k, v)| (*k, *v))
        .collect();

    for (key, range) in to_flush {
        flush_single_range(&range, key, node_iitree, path_iitree)?;
        range_buffer.remove(&key);
    }

    Ok(())
}

/// Break a large range into component ranges we haven't seen yet
///
/// Walk the range, breaking where we've seen it, emitting new ranges via lambda
fn for_each_fresh_range<F>(range: &Match, seen_bv: &[bool], mut lambda: F)
where
    F: FnMut(Match),
{
    let mut p = range.start;
    let mut t = range.data;

    while p < range.end {
        if seen_bv[p as usize] {
            p += 1;
            incr_pos(&mut t);
        } else {
            // Find the extent of the unseen range
            let q = p;
            let v = t;
            while p < range.end && !seen_bv[p as usize] {
                p += 1;
                incr_pos(&mut t);
            }
            lambda(Match::new(q, p, v));
        }
    }
}

/// Handle a range by marking it in the bitvector and adding to queues
fn handle_range(
    s: Match,
    curr_bv: &AtomicBitVec,
    ovlp_q: &OverlapAtomicQueue,
    todo_in: &RangeAtomicQueue,
) {
    let mut all_set_there = true;
    let mut n = s.data;
    for _i in s.start..s.end {
        // Set the bit and check if it was already set (using truly atomic fetch_or)
        let was_set = curr_bv.set(offset(n) as usize, true, Ordering::AcqRel);
        all_set_there = all_set_there && was_set;
        incr_pos(&mut n);
    }
    let _ = ovlp_q.push((s, is_rev(s.data)));
    if !all_set_there {
        let item = (make_pos_t(offset(s.data), is_rev(s.data)), s.end - s.start);
        let _ = todo_in.push(item);
    }
}

/// Explore overlaps from alignment iitree and handle them
fn explore_overlaps(
    b: &Match,
    seen_bv: &[bool],
    curr_bv: &AtomicBitVec,
    aln_iitree: &AdaptiveTree<u64, PosT>,
    ovlp_q: &OverlapAtomicQueue,
    todo_in: &RangeAtomicQueue,
) {
    aln_iitree
        .overlap(b.start, b.end, |_idx, start, end, pos| {
            let mut r = Match::new(start, end, pos);
            // Trim the range to fit within b
            if b.start > r.start {
                let trim_from_start = b.start - r.start;
                r.start += trim_from_start;
                incr_pos_by(&mut r.data, trim_from_start as usize);
            }
            if r.end > b.end {
                let trim_from_end = r.end - b.end;
                r.end -= trim_from_end;
            }
            assert!(r.start < r.end);
            for_each_fresh_range(&r, seen_bv, |s| {
                handle_range(s, curr_bv, ovlp_q, todo_in);
            });
        })
        .ok(); // Ignore Result
}

/// Write a chunk of the graph sequence from disjoint sets
fn write_graph_chunk(
    seqidx: &SeqIndex,
    node_iitree: &mut AdaptiveTree<u64, PosT>,
    path_iitree: &mut AdaptiveTree<u64, PosT>,
    seq_v_out: &mut Vec<u8>,
    range_buffer: &mut HashMap<PosT, Range>,
    dsets: Vec<(u64, u64)>,
    repeat_max: u64,
    min_repeat_dist: u64,
) -> io::Result<()> {
    let mut seq_v_length = seq_v_out.len() as u64;
    let mut last_dset_id = u64::MAX;
    let mut current_base = 0u8;

    let mut seq_counts: HashMap<u64, u64> = HashMap::new();
    let mut last_seq_pos: HashMap<u64, PosT> = HashMap::new();

    let close_to_prev = |seq_id: u64, pos: PosT, last_seq_pos: &HashMap<u64, PosT>| -> bool {
        if let Some(&last_pos) = last_seq_pos.get(&seq_id) {
            let dist = (offset(pos) as i64 - offset(last_pos) as i64).abs() as u64;
            dist < min_repeat_dist
        } else {
            false
        }
    };

    let mut todos: HashMap<u64, Vec<PosT>> = HashMap::new();

    for d in dsets {
        let curr_dset_id = d.0;
        let curr_offset = d.1;
        let base = seqidx.at(curr_offset).unwrap_or('N') as u8;

        // If we're on a new position
        if curr_dset_id != last_dset_id {
            if repeat_max != 0 || min_repeat_dist != 0 {
                // Flush todos inline
                for (_count, positions) in todos.iter() {
                    seq_v_out.push(current_base);
                    seq_v_length += 1;
                    for pos in positions {
                        extend_range(
                            seq_v_length - 1,
                            *pos,
                            range_buffer,
                            seqidx,
                            node_iitree,
                            path_iitree,
                        )?;
                    }
                }
                todos.clear();
                seq_counts.clear();
                last_seq_pos.clear();
            }
            // Emit new position
            current_base = base;
            seq_v_out.push(current_base);
            seq_v_length += 1;
            flush_ranges(seq_v_length - 1, range_buffer, node_iitree, path_iitree)?;
            last_dset_id = curr_dset_id;
        }

        let mut curr_q_pos = make_pos_t(curr_offset, false);
        if current_base != seqidx.at_pos(curr_q_pos).unwrap_or('N') as u8 {
            curr_q_pos = make_pos_t(curr_offset, true);
        }
        assert_eq!(current_base, seqidx.at_pos(curr_q_pos).unwrap_or('N') as u8);

        if let Some(curr_seq_id) = seqidx.seq_id_at(curr_offset) {
            let curr_seq_id = curr_seq_id as u64;
            let mut curr_seq_count = 0u64;

            if (min_repeat_dist != 0 && close_to_prev(curr_seq_id, curr_q_pos, &last_seq_pos))
                || (repeat_max != 0 && seq_counts.get(&curr_seq_id).unwrap_or(&0) + 1 > repeat_max)
            {
                curr_seq_count = *seq_counts.entry(curr_seq_id).or_insert(0) + 1;
                seq_counts.insert(curr_seq_id, curr_seq_count);
            } else if repeat_max != 0 || min_repeat_dist != 0 {
                *seq_counts.entry(curr_seq_id).or_insert(0) += 1;
            }

            if curr_seq_count == 0 {
                extend_range(
                    seq_v_length - 1,
                    curr_q_pos,
                    range_buffer,
                    seqidx,
                    node_iitree,
                    path_iitree,
                )?;
            } else {
                todos
                    .entry(curr_seq_count)
                    .or_insert_with(Vec::new)
                    .push(curr_q_pos);
            }
            last_seq_pos.insert(curr_seq_id, curr_q_pos);
        }
    }

    // Flush remaining todos
    for (_count, positions) in todos.iter() {
        seq_v_out.push(current_base);
        seq_v_length += 1;
        for pos in positions {
            extend_range(
                seq_v_length - 1,
                *pos,
                range_buffer,
                seqidx,
                node_iitree,
                path_iitree,
            )?;
        }
    }
    Ok(())
}

/// Main entry point for transitive closure computation
///
/// Computes connected components of aligned positions and builds
/// the variation graph sequence and interval trees.
pub fn compute_transitive_closures(
    seqidx: Arc<SeqIndex>,
    aln_iitree: Arc<AdaptiveTree<u64, PosT>>,
    seq_v_file: &str,
    node_iitree: Arc<RwLock<AdaptiveTree<u64, PosT>>>,
    path_iitree: Arc<RwLock<AdaptiveTree<u64, PosT>>>,
    repeat_max: u64,
    min_repeat_dist: u64,
    transclose_batch_size: u64,
    show_progress: bool,
    num_threads: usize,
) -> io::Result<usize> {
    use std::collections::VecDeque;
    use std::fs::File;
    use std::io::Write;

    let start_time = std::time::Instant::now();
    eprintln!("[transclosure] Starting transitive closure computation");
    eprintln!("[transclosure] Using {} threads", num_threads);

    // Open iitree writers (need write access)
    node_iitree.write().unwrap().open_writer()?;
    path_iitree.write().unwrap().open_writer()?;

    // Open output file
    let mut seq_v_out = Some(Vec::new());

    // Bitvector to track visited positions
    let input_seq_length = seqidx.seq_length() as usize;
    let mut q_seen_bv = vec![false; input_seq_length];

    // Range buffer for writing to iitrees
    let mut range_buffer = Some(HashMap::new());

    let mut bases_seen = 0u64;

    // Writer thread handle for pipelining (like C++)
    // The thread writes while we compute the next batch
    let mut writer_thread: Option<thread::JoinHandle<io::Result<(Vec<u8>, HashMap<PosT, Range>)>>> =
        None;

    // Main loop: process input sequence in chunks
    let mut i = 0;
    while i < input_seq_length {
        // Skip already-seen positions
        while i < input_seq_length && q_seen_bv[i] {
            i += 1;
        }
        if i >= input_seq_length {
            break;
        }

        let chunk_start = i;
        let mut bases_to_consider = 0;
        let mut chunk_end = chunk_start;
        while bases_to_consider < transclose_batch_size as usize && chunk_end < input_seq_length {
            if !q_seen_bv[chunk_end] {
                bases_to_consider += 1;
            }
            chunk_end += 1;
        }

        if show_progress {
            eprintln!(
                "[transclosure] {:.3}s {:.2}% {}-{} overlap_collect",
                start_time.elapsed().as_secs_f64(),
                (bases_seen as f64 / input_seq_length as f64) * 100.0,
                chunk_start,
                chunk_end
            );
        }

        // Atomic bitvector for positions seen in this chunk
        let q_curr_bv = AtomicBitVec::new(input_seq_length);

        // Work queues
        let todo_in = Arc::new(ArrayQueue::new(100000));
        let todo_out = Arc::new(ArrayQueue::new(100000));
        let ovlp_q = Arc::new(ArrayQueue::new(100000));

        let mut todo: VecDeque<(PosT, u64)> = VecDeque::new();
        let mut ovlp: Vec<(Match, bool)> = Vec::new();

        // Seed initial ranges
        for_each_fresh_range(
            &Match::new(chunk_start as u64, chunk_end as u64, 0),
            &q_seen_bv,
            |b| {
                for j in b.start..b.end {
                    q_curr_bv.set(j as usize, true, Ordering::Release);
                }
                let range = (make_pos_t(b.start, false), b.end - b.start);
                if todo_out.push(range).is_err() {
                    todo.push_back(range);
                }
            },
        );

        // Parallel exploration using Rayon's work-stealing
        // This eliminates the manager thread bottleneck and thread oversubscription
        let aln_iitree_clone = Arc::clone(&aln_iitree);
        let q_curr_bv_shared = Arc::new(q_curr_bv);

        // Use Rayon's scope for dynamic task spawning with work-stealing
        rayon::scope(|s| {
            // Spawn worker tasks (2x threads for better load balancing)
            for _ in 0..(num_threads * 2) {
                let todo_out = Arc::clone(&todo_out);
                let todo_in = Arc::clone(&todo_in);
                let ovlp_q = Arc::clone(&ovlp_q);
                let aln_iitree = Arc::clone(&aln_iitree_clone);
                let q_curr_bv = Arc::clone(&q_curr_bv_shared);
                let q_seen_bv_clone = q_seen_bv.clone();

                s.spawn(move |_s| {
                    // Process work items until queue is empty
                    loop {
                        // Try to get work from todo_out
                        if let Some((pos, match_len)) = todo_out.pop() {
                            let n = if !is_rev(pos) {
                                offset(pos)
                            } else {
                                offset(pos) - match_len + 1
                            };
                            let range_start = n;
                            let range_end = n + match_len;

                            explore_overlaps(
                                &Match::new(range_start, range_end, pos),
                                &q_seen_bv_clone,
                                &q_curr_bv,
                                &aln_iitree,
                                &ovlp_q,
                                &todo_in,
                            );
                        } else {
                            // No work available, yield to scheduler (like C++ 1ns sleep)
                            std::thread::yield_now();

                            // If still no work and no work in todo_in, we're likely done
                            if todo_out.is_empty() && todo_in.is_empty() {
                                break;
                            }
                        }
                    }
                });
            }

            // Manager task to move work from todo_in -> todo -> todo_out
            // This runs in parallel with workers on Rayon's thread pool
            s.spawn(move |_s| {
                let mut empty_count = 0;
                loop {
                    let mut did_work = false;

                    // Transfer from todo_in to todo
                    while let Some(item) = todo_in.pop() {
                        todo.push_back(item);
                        did_work = true;
                    }

                    // Transfer from todo to todo_out
                    while let Some(item) = todo.front().copied() {
                        if todo_out.push(item).is_ok() {
                            todo.pop_front();
                            did_work = true;
                        } else {
                            break;
                        }
                    }

                    if did_work {
                        empty_count = 0;
                    } else {
                        std::thread::yield_now();
                        empty_count += 1;
                        // If queues have been empty for a while, we're done
                        if empty_count > 1000
                            && todo.is_empty()
                            && todo_in.is_empty()
                            && todo_out.is_empty()
                        {
                            break;
                        }
                    }
                }
            });
        });

        // Collect all overlaps after workers complete
        while let Some(o) = ovlp_q.pop() {
            ovlp.push(o);
        }

        if show_progress {
            eprintln!(
                "[transclosure] {:.3}s {:.2}% {}-{} union_find",
                start_time.elapsed().as_secs_f64(),
                (bases_seen as f64 / input_seq_length as f64) * 100.0,
                chunk_start,
                chunk_end
            );
        }

        // Build dense mapping using rank (like C++ sdsl::rank_1_type)
        if show_progress {
            eprintln!(
                "[transclosure] {:.3}s {:.2}% {}-{} rank_build_start",
                start_time.elapsed().as_secs_f64(),
                (bases_seen as f64 / input_seq_length as f64) * 100.0,
                chunk_start,
                chunk_end
            );
        }

        let q_curr_bv_final = Arc::try_unwrap(q_curr_bv_shared).unwrap();

        // Parallelize rank building (like C++ parallel_for for q_curr_bv_vec)
        // Process in parallel chunks of 100000 positions
        // Collect only the set bit positions in parallel (no need for full Vec<bool>)
        let chunk_size = 100000;
        let num_chunks = (input_seq_length + chunk_size - 1) / chunk_size;

        let q_curr_positions: Vec<u64> = (0..num_chunks)
            .into_par_iter()
            .flat_map(|chunk_idx| {
                let start = chunk_idx * chunk_size;
                let end = (start + chunk_size).min(input_seq_length);
                let mut local_positions = Vec::new();

                for pos in start..end {
                    if q_curr_bv_final.get(pos, Ordering::Acquire) {
                        local_positions.push(pos as u64);
                    }
                }
                local_positions
            })
            .collect();

        let q_curr_bv_count = q_curr_positions.len();
        if q_curr_bv_count == 0 {
            i = chunk_end;
            continue;
        }

        // Build Rank9Sel directly from iterator (avoids intermediate Vec<bool> allocation)
        let sucds_bv = SucdsBitVector::from_bits(
            (0..input_seq_length).map(|pos| q_curr_bv_final.get(pos, Ordering::Acquire)),
        );
        let q_curr_rank = Rank9Sel::new(sucds_bv);

        if show_progress {
            eprintln!(
                "[transclosure] {:.3}s {:.2}% {}-{} rank_build_end",
                start_time.elapsed().as_secs_f64(),
                (bases_seen as f64 / input_seq_length as f64) * 100.0,
                chunk_start,
                chunk_end
            );
        }

        // Union-find with O(1) rank lookups (using INLINE ASM to match C++ perf)
        let dsets = DisjointSetsAsm::new(q_curr_bv_count);

        // Use chunked iteration with grain size 10000 (like C++ paryfor grain size)
        // This reduces contention on the lock-free union-find structure
        ovlp.par_chunks(10000).for_each(|chunk| {
            for s in chunk {
                let r = &s.0;
                let mut p = r.data;
                for j in r.start..r.end {
                    // O(1) rank lookups: rank1(pos) counts 1-bits in [0..pos), giving 0-indexed rank
                    let j_rank = q_curr_rank.rank1(j as usize).unwrap();
                    let p_rank = q_curr_rank.rank1(offset(p) as usize).unwrap();
                    dsets.unite(j_rank, p_rank);
                    incr_pos(&mut p);
                }
            }
        });

        if show_progress {
            eprintln!(
                "[transclosure] {:.3}s {:.2}% {}-{} dset_write",
                start_time.elapsed().as_secs_f64(),
                (bases_seen as f64 / input_seq_length as f64) * 100.0,
                chunk_start,
                chunk_end
            );
        }

        // Read out disjoint sets
        // Use chunked iteration with grain size 10000 to reduce overhead
        let mut dsets_vec: Vec<(u64, u64)> = (0..q_curr_positions.len())
            .into_par_iter()
            .with_min_len(10000)
            .filter_map(|j| {
                let p = q_curr_positions[j];
                if !q_seen_bv[p as usize] {
                    Some((dsets.find(j) as u64, p))
                } else {
                    None
                }
            })
            .collect();

        if dsets_vec.is_empty() {
            i = chunk_end;
            continue;
        }

        // Sort and compress
        if show_progress {
            eprintln!(
                "[transclosure] {:.3}s {:.2}% {}-{} dset_sort_start",
                start_time.elapsed().as_secs_f64(),
                (bases_seen as f64 / input_seq_length as f64) * 100.0,
                chunk_start,
                chunk_end
            );
        }

        dsets_vec.par_sort_unstable();

        if show_progress {
            eprintln!(
                "[transclosure] {:.3}s {:.2}% {}-{} dset_sort1_end",
                start_time.elapsed().as_secs_f64(),
                (bases_seen as f64 / input_seq_length as f64) * 100.0,
                chunk_start,
                chunk_end
            );
        }

        // Compress dset IDs
        let mut c = 0u64;
        let mut last_id = dsets_vec[0].0;
        for d in &mut dsets_vec {
            if d.0 != last_id {
                c += 1;
                last_id = d.0;
            }
            d.0 = c;
        }

        // Find minimum position in each dset
        let mut dsets_by_min_pos = vec![(u64::MAX, 0u64); (c + 1) as usize];
        for i in 0..=c {
            dsets_by_min_pos[i as usize].1 = i;
        }
        for d in &dsets_vec {
            let minpos = &mut dsets_by_min_pos[d.0 as usize].0;
            *minpos = (*minpos).min(d.1);
        }

        if show_progress {
            eprintln!(
                "[transclosure] {:.3}s {:.2}% {}-{} dset_sort2_start",
                start_time.elapsed().as_secs_f64(),
                (bases_seen as f64 / input_seq_length as f64) * 100.0,
                chunk_start,
                chunk_end
            );
        }

        dsets_by_min_pos.par_sort_unstable();

        if show_progress {
            eprintln!(
                "[transclosure] {:.3}s {:.2}% {}-{} dset_sort2_end",
                start_time.elapsed().as_secs_f64(),
                (bases_seen as f64 / input_seq_length as f64) * 100.0,
                chunk_start,
                chunk_end
            );
        }

        // Invert naming
        let mut dset_names = vec![0u64; (c + 1) as usize];
        for (x, d) in dsets_by_min_pos.iter().enumerate() {
            dset_names[d.1 as usize] = x as u64;
        }

        // Rename and re-sort
        for d in &mut dsets_vec {
            d.0 = dset_names[d.0 as usize];
        }

        if show_progress {
            eprintln!(
                "[transclosure] {:.3}s {:.2}% {}-{} dset_sort3_start",
                start_time.elapsed().as_secs_f64(),
                (bases_seen as f64 / input_seq_length as f64) * 100.0,
                chunk_start,
                chunk_end
            );
        }

        dsets_vec.par_sort_unstable();

        if show_progress {
            eprintln!(
                "[transclosure] {:.3}s {:.2}% {}-{} dset_sort3_end",
                start_time.elapsed().as_secs_f64(),
                (bases_seen as f64 / input_seq_length as f64) * 100.0,
                chunk_start,
                chunk_end
            );
        }

        // Mark as seen
        for d in &dsets_vec {
            q_seen_bv[d.1 as usize] = true;
            bases_seen += 1;
        }

        if show_progress {
            eprintln!(
                "[transclosure] {:.3}s {:.2}% {}-{} graph_emission",
                start_time.elapsed().as_secs_f64(),
                (bases_seen as f64 / input_seq_length as f64) * 100.0,
                chunk_start,
                chunk_end
            );
        }

        // Pipeline write_graph_chunk (like C++):
        // Join previous writer (if any) to get back seq_v_out and range_buffer
        if let Some(handle) = writer_thread.take() {
            let result = handle.join().expect("Writer thread panicked")?;
            seq_v_out = Some(result.0);
            range_buffer = Some(result.1);
        }

        // Take ownership of seq_v_out and range_buffer for the thread
        let mut seq_v = seq_v_out.take().unwrap();
        let mut range_buf = range_buffer.take().unwrap();

        // Spawn new writer thread with cloned Arc handles
        let node_iitree_clone = Arc::clone(&node_iitree);
        let path_iitree_clone = Arc::clone(&path_iitree);
        let seqidx_clone = Arc::clone(&seqidx);

        writer_thread = Some(thread::spawn(
            move || -> io::Result<(Vec<u8>, HashMap<PosT, Range>)> {
                let mut node_guard = node_iitree_clone.write().unwrap();
                let mut path_guard = path_iitree_clone.write().unwrap();

                write_graph_chunk(
                    &seqidx_clone,
                    &mut node_guard,
                    &mut path_guard,
                    &mut seq_v,
                    &mut range_buf,
                    dsets_vec,
                    repeat_max,
                    min_repeat_dist,
                )?;

                Ok((seq_v, range_buf))
            },
        ));

        i = chunk_end;
    }

    // Join the final writer thread
    if let Some(handle) = writer_thread.take() {
        let result = handle.join().expect("Writer thread panicked")?;
        seq_v_out = Some(result.0);
        range_buffer = Some(result.1);
    }

    // Unwrap the final values
    let seq_v_out = seq_v_out.unwrap();
    let mut range_buffer = range_buffer.unwrap();

    // Write output file
    let mut file = File::create(seq_v_file)?;
    file.write_all(&seq_v_out)?;
    let seq_bytes = seq_v_out.len();

    // Flush remaining ranges
    {
        let mut node_guard = node_iitree.write().unwrap();
        let mut path_guard = path_iitree.write().unwrap();
        flush_ranges(
            seq_bytes as u64 + 1,
            &mut range_buffer,
            &mut node_guard,
            &mut path_guard,
        )?;
    }

    if show_progress {
        eprintln!("[transclosure] Building node_iitree and path_iitree indexes");
    }

    // Close writers and build indexes
    node_iitree.write().unwrap().close_writer()?;
    path_iitree.write().unwrap().close_writer()?;
    node_iitree.write().unwrap().index()?;
    path_iitree.write().unwrap().index()?;

    eprintln!("[transclosure] Transitive closure computation complete");

    Ok(seq_bytes)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_wang_hash() {
        // Wang hash should be deterministic
        let val = 12345u64;
        assert_eq!(wang_hash_64(val), wang_hash_64(val));

        // Different values should hash differently
        assert_ne!(wang_hash_64(12345), wang_hash_64(54321));
    }

    #[test]
    fn test_range() {
        let r = Range::new(10, 20);
        assert_eq!(r.begin, 10);
        assert_eq!(r.end, 20);
    }

    #[test]
    fn test_match() {
        let m = Match::new(0, 100, make_pos_t(50, false));
        assert_eq!(m.length(), 100);
    }
}
