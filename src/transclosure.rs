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

use crate::intervaltree::AdaptiveTree;
use crate::intervaltree::IntervalTree;
use crate::pos::{decr_pos, decr_pos_by, incr_pos, incr_pos_by, is_rev, make_pos_t, offset, PosT};
use crate::seqindex::SeqIndex;

const BFS_QUEUE_CAPACITY: usize = 1 << 17;
const RANK_CHUNK_SIZE: usize = 1 << 17;

enum UnionFindResult {
    Gpu(Vec<u32>),
    Cpu(DisjointSetsAsm),
}

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

/// Type alias for BFS work queue
type RangeAtomicQueue = ArrayQueue<(PosT, u64)>;

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
    node_iitree.add(match_start_in_s, match_end_in_s, match_pos_in_q)?;
    path_iitree.add(match_start_in_q, match_end_in_q, match_pos_in_s)?;

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

/// Find sequences with at least one discovered position in the bitvector.
fn find_component_sequences(seqidx: &SeqIndex, bv: &AtomicBitVec) -> Vec<usize> {
    (1..=seqidx.n_seqs())
        .filter(|&seq_id| {
            if let Some(off) = seqidx.nth_seq_offset(seq_id) {
                bv.get(off as usize, Ordering::Relaxed)
            } else {
                false
            }
        })
        .collect()
}

/// Explore overlaps from alignment iitree — discovery only, no ovlp_q collection.
/// Only follows edges in the spanning tree for fast BFS.
fn explore_overlaps_discovery(
    b: &Match,
    seen_bv: &[bool],
    curr_bv: &AtomicBitVec,
    aln_iitree: &AdaptiveTree<u64, PosT>,
    todo_in: &RangeAtomicQueue,
    seqidx: &SeqIndex,
    spanning_adj: &SpanningTreeAdj,
) {
    let source_seq = seqidx.seq_id_at(b.start).unwrap_or(0);

    aln_iitree
        .overlap(b.start, b.end, |_idx, start, end, pos| {
            // Filter: only follow spanning tree edges
            let target_seq = seqidx.seq_id_at(offset(pos)).unwrap_or(0);
            if !spanning_adj.contains(source_seq, target_seq) {
                return;
            }

            let mut r = Match::new(start, end, pos);
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
                // Discovery only: set curr_bv bits, push to todo_in if new
                // No ovlp_q push — union-find is handled in phase 2
                let mut all_set_there = true;
                let mut n = s.data;
                for _i in s.start..s.end {
                    let was_set = curr_bv.set(offset(n) as usize, true, Ordering::AcqRel);
                    all_set_there = all_set_there && was_set;
                    incr_pos(&mut n);
                }
                if !all_set_there {
                    let item = (make_pos_t(offset(s.data), is_rev(s.data)), s.end - s.start);
                    while todo_in.push(item).is_err() {
                        std::thread::yield_now();
                    }
                }
            });
        })
        .ok();
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
            let dist = (offset(pos) as i64 - offset(last_pos) as i64).unsigned_abs();
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
                todos.entry(curr_seq_count).or_default().push(curr_q_pos);
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

/// Adjacency list for spanning tree: adj[seq_id] contains the neighbor seq_ids.
/// Supports O(degree) lookup — much faster than HashSet for small degree (~2-3).
struct SpanningTreeAdj {
    adj: Vec<Vec<usize>>,
}

impl SpanningTreeAdj {
    fn new(n_seqs: usize) -> Self {
        SpanningTreeAdj {
            adj: vec![Vec::new(); n_seqs + 1],
        }
    }

    fn add_edge(&mut self, a: usize, b: usize) {
        self.adj[a].push(b);
        self.adj[b].push(a);
    }

    #[inline]
    fn contains(&self, source: usize, target: usize) -> bool {
        // Degree is ~2-3 for a tree, so linear scan is faster than hashing
        self.adj[source].contains(&target)
    }
}

/// Compute a maximum-weight spanning tree of sequence pairs from the alignment iitree.
///
/// Returns `(adj, tree_edges, total_pairs)` for use in progress logging.
fn compute_spanning_tree(
    aln_iitree: &AdaptiveTree<u64, PosT>,
    seqidx: &SeqIndex,
) -> (SpanningTreeAdj, usize, usize) {
    let n_seqs = seqidx.n_seqs();
    if n_seqs <= 1 {
        return (SpanningTreeAdj::new(n_seqs), 0, 0);
    }

    // Collect pair weights: total aligned bases per (seq_i, seq_j) pair
    let mut pair_weights: HashMap<(usize, usize), u64> = HashMap::new();
    aln_iitree
        .for_each_interval(|start, end, target_pos| {
            let source_seq = seqidx.seq_id_at(start).unwrap_or(0);
            let target_seq = seqidx.seq_id_at(offset(target_pos)).unwrap_or(0);
            if source_seq != target_seq && source_seq > 0 && target_seq > 0 {
                let canonical = if source_seq <= target_seq {
                    (source_seq, target_seq)
                } else {
                    (target_seq, source_seq)
                };
                *pair_weights.entry(canonical).or_insert(0) += end - start;
            }
        })
        .ok();

    // Sort pairs by weight descending (max-weight spanning tree)
    let mut edges: Vec<((usize, usize), u64)> = pair_weights.into_iter().collect();
    edges.sort_unstable_by_key(|edge| std::cmp::Reverse(edge.1));

    // Kruskal's algorithm with simple union-find
    let mut parent: Vec<usize> = (0..=n_seqs).collect();
    let mut rank: Vec<usize> = vec![0; n_seqs + 1];

    fn find(parent: &mut [usize], x: usize) -> usize {
        if parent[x] != x {
            parent[x] = find(parent, parent[x]);
        }
        parent[x]
    }

    fn unite(parent: &mut [usize], rank: &mut [usize], x: usize, y: usize) -> bool {
        let rx = find(parent, x);
        let ry = find(parent, y);
        if rx == ry {
            return false;
        }
        if rank[rx] < rank[ry] {
            parent[rx] = ry;
        } else if rank[rx] > rank[ry] {
            parent[ry] = rx;
        } else {
            parent[ry] = rx;
            rank[rx] += 1;
        }
        true
    }

    let mut spanning_adj = SpanningTreeAdj::new(n_seqs);
    let mut tree_edges = 0;

    for ((s1, s2), _weight) in &edges {
        if unite(&mut parent, &mut rank, *s1, *s2) {
            spanning_adj.add_edge(*s1, *s2);
            tree_edges += 1;
            if tree_edges >= n_seqs - 1 {
                break;
            }
        }
    }

    (spanning_adj, tree_edges, edges.len())
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
    if show_progress {
        eprintln!("[transclosure] Starting transitive closure computation");
        eprintln!("[transclosure] Using {num_threads} threads");
    }

    // Compute spanning tree for fast BFS discovery
    let (spanning_pairs, st_edges, st_pairs) = compute_spanning_tree(&aln_iitree, &seqidx);
    if show_progress {
        eprintln!(
            "[transclosure] Spanning tree: {st_edges} edges from {st_pairs} total pairs ({}x reduction)",
            st_pairs.checked_div(st_edges).unwrap_or(0)
        );
    }

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

    // Initialise the GPU runner once here. Reused across all batches so
    // the CUDA context, PTX module, and pinned/device buffers are just allocated once. Falls back to CPU if CUDA is not available.
    let mut gpu_runner: Option<crate::gpu::GpuRunner> = crate::gpu::GpuRunner::new();

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

        // ================================================================
        // PHASE 1: Fast BFS discovery using spanning tree edges only
        // ================================================================
        // Only follows spanning tree edges (N-1 pairs instead of N*(N-1)/2).
        // No overlap collection (ovlp_q) — just discovers positions in curr_bv.
        let q_curr_bv = AtomicBitVec::new(input_seq_length);

        // Work queues for BFS (no ovlp_q needed)
        let todo_in = Arc::new(ArrayQueue::new(BFS_QUEUE_CAPACITY));
        let todo_out = Arc::new(ArrayQueue::new(BFS_QUEUE_CAPACITY));
        let mut todo: VecDeque<(PosT, u64)> = VecDeque::new();
        let active_workers = Arc::new(AtomicU64::new(0));

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

        // Parallel BFS — discovery only, spanning tree edges only
        let aln_iitree_clone = Arc::clone(&aln_iitree);
        let seqidx_clone = Arc::clone(&seqidx);
        let q_curr_bv_shared = Arc::new(q_curr_bv);
        let spanning_pairs_ref = &spanning_pairs;

        rayon::scope(|s| {
            for _worker_idx in 0..(num_threads * 2) {
                let todo_out = Arc::clone(&todo_out);
                let todo_in = Arc::clone(&todo_in);
                let aln_iitree = Arc::clone(&aln_iitree_clone);
                let q_curr_bv = Arc::clone(&q_curr_bv_shared);
                let q_seen_bv_clone = q_seen_bv.clone();
                let seqidx = Arc::clone(&seqidx_clone);

                let active_workers_clone = Arc::clone(&active_workers);
                s.spawn(move |_s| {
                    let mut empty_count = 0u64;
                    loop {
                        if let Some((pos, match_len)) = todo_out.pop() {
                            empty_count = 0;
                            active_workers_clone.fetch_add(1, Ordering::SeqCst);
                            let n = if !is_rev(pos) {
                                offset(pos)
                            } else {
                                offset(pos) - match_len + 1
                            };

                            explore_overlaps_discovery(
                                &Match::new(n, n + match_len, pos),
                                &q_seen_bv_clone,
                                &q_curr_bv,
                                &aln_iitree,
                                &todo_in,
                                &seqidx,
                                spanning_pairs_ref,
                            );
                            active_workers_clone.fetch_sub(1, Ordering::SeqCst);
                        } else {
                            std::thread::yield_now();
                            empty_count += 1;
                            let to_empty = todo_out.is_empty();
                            let ti_empty = todo_in.is_empty();
                            let no_active = active_workers_clone.load(Ordering::SeqCst) == 0;
                            if to_empty && ti_empty && no_active {
                                if empty_count > 1000 {
                                    break;
                                }
                            } else if !no_active {
                                empty_count = 0;
                            }
                        }
                    }
                });
            }

            // Manager task — shuttles todo_in → todo → todo_out (no ovlp_q to drain)
            let active_workers_mgr = Arc::clone(&active_workers);
            s.spawn(move |_s| {
                let mut empty_count = 0;
                loop {
                    let mut did_work = false;
                    while let Some(item) = todo_in.pop() {
                        todo.push_back(item);
                        did_work = true;
                    }
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
                        let no_active = active_workers_mgr.load(Ordering::SeqCst) == 0;
                        if empty_count > 1000
                            && todo.is_empty()
                            && todo_in.is_empty()
                            && todo_out.is_empty()
                            && no_active
                        {
                            break;
                        }
                        if !no_active {
                            empty_count = 0;
                        }
                    }
                }
            });
        });

        if show_progress {
            eprintln!(
                "[transclosure] {:.3}s {:.2}% {}-{} phase1_discovery_done",
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

        // ================================================================
        // PHASE 1b: Orphan recovery — flood-fill using per-sequence queries
        // ================================================================
        // The spanning tree BFS may miss positions only reachable via non-tree
        // edges. Query the iitree per-sequence (only sequences with discovered
        // positions) to find and mark orphans. Converges when no new positions found.
        {
            let mut recovery_round = 0u32;
            loop {
                let new_positions = AtomicU64::new(0);
                // Find sequences with discovered positions
                let component_seqs = find_component_sequences(&seqidx, &q_curr_bv_final);

                component_seqs.par_iter().for_each(|&seq_id| {
                    let seq_off = seqidx.nth_seq_offset(seq_id).unwrap();
                    let seq_len = seqidx.nth_seq_length(seq_id).unwrap();
                    aln_iitree
                        .overlap(
                            seq_off,
                            seq_off + seq_len,
                            |_idx, iv_start, iv_end, target_pos| {
                                if !q_curr_bv_final.get(iv_start as usize, Ordering::Relaxed) {
                                    return;
                                }
                                let mut p = target_pos;
                                for _ in iv_start..iv_end {
                                    let t = offset(p) as usize;
                                    if t < input_seq_length
                                        && !q_curr_bv_final.get(t, Ordering::Relaxed)
                                    {
                                        let was_set =
                                            q_curr_bv_final.set(t, true, Ordering::AcqRel);
                                        if !was_set {
                                            new_positions.fetch_add(1, Ordering::Relaxed);
                                        }
                                    }
                                    incr_pos(&mut p);
                                }
                            },
                        )
                        .ok();
                });
                let found = new_positions.load(Ordering::Relaxed);
                if found == 0 {
                    break;
                }
                recovery_round += 1;
                eprintln!(
                    "[transclosure] {:.3}s orphan_recovery round {}: {} new positions ({} seqs in component)",
                    start_time.elapsed().as_secs_f64(),
                    recovery_round,
                    found,
                    component_seqs.len(),
                );
            }
            if recovery_round > 0 {
                eprintln!(
                    "[transclosure] {:.3}s orphan_recovery complete after {} rounds",
                    start_time.elapsed().as_secs_f64(),
                    recovery_round
                );
            }
        }

        // Collect set-bit positions in sorted order (flat_map preserves chunk order).
        let num_chunks = (input_seq_length + RANK_CHUNK_SIZE - 1) / RANK_CHUNK_SIZE;

        let q_curr_positions: Vec<u64> = (0..num_chunks)
            .into_par_iter()
            .flat_map(|chunk_idx| {
                let start = chunk_idx * RANK_CHUNK_SIZE;
                let end = (start + RANK_CHUNK_SIZE).min(input_seq_length);
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

        // Position-to-rank scatter. q_curr_positions is sorted, so writes are
        // sequential and cache-friendly.
        let mut rank_table = vec![0u32; input_seq_length];
        for (rank, &pos) in q_curr_positions.iter().enumerate() {
            rank_table[pos as usize] = rank as u32;
        }

        if show_progress {
            eprintln!(
                "[transclosure] {:.3}s {:.2}% {}-{} rank_build_end",
                start_time.elapsed().as_secs_f64(),
                (bases_seen as f64 / input_seq_length as f64) * 100.0,
                chunk_start,
                chunk_end
            );
        }

        let q_curr_bv_ref = &q_curr_bv_final;

        if show_progress {
            eprintln!(
                "[transclosure] {:.3}s {:.2}% {}-{} phase2_union_find_start",
                start_time.elapsed().as_secs_f64(),
                (bases_seen as f64 / input_seq_length as f64) * 100.0,
                chunk_start,
                chunk_end
            );
        }

        let component_seqs = find_component_sequences(&seqidx, &q_curr_bv_final);

        // CUDA availability is stable for the process lifetime; check once.
        use std::sync::OnceLock;
        static CUDA_AVAILABLE: OnceLock<bool> = OnceLock::new();
        // Dynamic crossover threshold (3M elements) based on benchmarks where GPU outperforms CPU
        const GPU_MIN_ELEMENTS_THRESHOLD: usize = 3_000_000;
        let force_gpu = std::env::var("SEQWISH_FORCE_GPU").is_ok();
        let use_gpu = *CUDA_AVAILABLE.get_or_init(crate::gpu::is_cuda_available)
            && (force_gpu || q_curr_bv_count >= GPU_MIN_ELEMENTS_THRESHOLD);

        // `gpu_runner` is initialised once before the main loop (below) and
        // reused across batches so the CUDA context, PTX module, and device
        // buffers are not re-created on every chunk.

        let uf_result: UnionFindResult = if use_gpu {
            let edges: Vec<(u32, u32)> = component_seqs
                .par_iter()
                .flat_map(|&seq_id| {
                    let seq_off = seqidx.nth_seq_offset(seq_id).unwrap();
                    let seq_len = seqidx.nth_seq_length(seq_id).unwrap();
                    let mut local_edges = Vec::new();
                    let _ = aln_iitree.overlap(
                        seq_off,
                        seq_off + seq_len,
                        |_idx, start, end, target_pos| {
                            if !q_curr_bv_ref.get(start as usize, Ordering::Relaxed) {
                                return;
                            }
                            let mut p = target_pos;
                            for j in start..end {
                                let t = offset(p) as usize;
                                if t < input_seq_length
                                    && q_curr_bv_ref.get(j as usize, Ordering::Relaxed)
                                    && q_curr_bv_ref.get(t, Ordering::Relaxed)
                                {
                                    local_edges.push((rank_table[j as usize], rank_table[t]));
                                }
                                incr_pos(&mut p);
                            }
                        },
                    );
                    local_edges
                })
                .collect();

            match gpu_runner
                .as_mut()
                .and_then(|r| r.gpu_union_find(q_curr_bv_count, &edges, show_progress))
            {
                Some(roots) => UnionFindResult::Gpu(roots),
                None => {
                    // check succeeded but context creation failed (no device or OOM).
                    let cpu = DisjointSetsAsm::new(q_curr_bv_count);
                    edges.par_iter().for_each(|&(u, v)| {
                        cpu.unite(u as usize, v as usize);
                    });
                    UnionFindResult::Cpu(cpu)
                }
            }
        } else {
            let cpu = DisjointSetsAsm::new(q_curr_bv_count);
            component_seqs.par_iter().for_each(|&seq_id| {
                let seq_off = seqidx.nth_seq_offset(seq_id).unwrap();
                let seq_len = seqidx.nth_seq_length(seq_id).unwrap();
                let _ = aln_iitree.overlap(
                    seq_off,
                    seq_off + seq_len,
                    |_idx, start, end, target_pos| {
                        if !q_curr_bv_ref.get(start as usize, Ordering::Relaxed) {
                            return;
                        }
                        let mut p = target_pos;
                        for j in start..end {
                            let t = offset(p) as usize;
                            if t < input_seq_length
                                && q_curr_bv_ref.get(j as usize, Ordering::Relaxed)
                                && q_curr_bv_ref.get(t, Ordering::Relaxed)
                            {
                                cpu.unite(rank_table[j as usize] as usize, rank_table[t] as usize);
                            }
                            incr_pos(&mut p);
                        }
                    },
                );
            });
            UnionFindResult::Cpu(cpu)
        };

        if show_progress {
            eprintln!(
                "[transclosure] {:.3}s {:.2}% {}-{} phase2_union_find_done ({} seqs)",
                start_time.elapsed().as_secs_f64(),
                (bases_seen as f64 / input_seq_length as f64) * 100.0,
                chunk_start,
                chunk_end,
                component_seqs.len(),
            );
        }

        if show_progress {
            eprintln!(
                "[transclosure] {:.3}s {:.2}% {}-{} dset_write",
                start_time.elapsed().as_secs_f64(),
                (bases_seen as f64 / input_seq_length as f64) * 100.0,
                chunk_start,
                chunk_end
            );
        }

        // Map each rank to its root, along with its input position.
        let mut dsets_vec: Vec<(u64, u64)> = (0..q_curr_positions.len())
            .into_par_iter()
            .with_min_len(10_000)
            .filter_map(|j| {
                let p = q_curr_positions[j];
                if q_seen_bv[p as usize] {
                    return None;
                }
                let root = match &uf_result {
                    UnionFindResult::Gpu(roots) => roots[j] as u64,
                    UnionFindResult::Cpu(cpu) => cpu.find(j) as u64,
                };
                Some((root, p))
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

    // Close writers and build indexes
    if show_progress {
        eprintln!("[transclosure] Building node_iitree and path_iitree indexes");
    }
    node_iitree.write().unwrap().close_writer()?;
    path_iitree.write().unwrap().close_writer()?;
    node_iitree.write().unwrap().index()?;
    path_iitree.write().unwrap().index()?;

    if show_progress {
        eprintln!("[transclosure] Transitive closure computation complete");
    }

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
