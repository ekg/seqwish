#pragma once

#include <string>
#include <fstream>
#include <iostream>
#include <unordered_set>
#include <set>
#include <thread>
#include "sdsl/bit_vectors.hpp"
#include "atomic_bitvector.hpp"
#include "flat_hash_map.hpp"
#include "seqindex.hpp"
#include "mmiitree.hpp"
#include "pos.hpp"
#include "match.hpp"
#include "ips4o.hpp"
#include "spinlock.hpp"
#include "dset64-gccAtomic.hpp"
#include "atomic_queue.h"
#include "time.hpp"
#include "wang.hpp"

namespace seqwish {

#define DEBUG_TRANSCLOSURE true

typedef atomic_queue::AtomicQueue2<std::pair<pos_t, uint64_t>, 2 << 16> range_atomic_queue_t;
typedef atomic_queue::AtomicQueue2<std::pair<match_t, bool>, 2 << 16> overlap_atomic_queue_t;

struct range_t {
    uint64_t begin = 0;
    uint64_t end = 0;
};

void extend_range(const uint64_t& s_pos,
                  const pos_t& q_pos,
                  std::map<pos_t, range_t>& range_buffer,
                  const seqindex_t& seqidx,
                  mmmulti::iitree<uint64_t, pos_t>& node_iitree,
                  mmmulti::iitree<uint64_t, pos_t>& path_iitree);

void flush_ranges(const uint64_t& s_pos,
                  std::map<pos_t, range_t>& range_buffer,
                  mmmulti::iitree<uint64_t, pos_t>& node_iitree,
                  mmmulti::iitree<uint64_t, pos_t>& path_iitree);

void flush_range(std::map<pos_t, range_t>::iterator it,
                 mmmulti::iitree<uint64_t, pos_t>& node_iitree,
                 mmmulti::iitree<uint64_t, pos_t>& path_iitree);

void for_each_fresh_range(const match_t& range,
                          const std::vector<bool>& seen_bv,
                          const std::function<void(match_t)>& lambda);

void handle_range(match_t s,
                  atomicbitvector::atomic_bv_t& curr_bv,
                  overlap_atomic_queue_t& ovlp_q,
                  range_atomic_queue_t& todo_in);

void explore_overlaps(const match_t& b,
                      const std::vector<bool>& seen_bv,
                      atomicbitvector::atomic_bv_t& curr_bv,
                      mmmulti::iitree<uint64_t, pos_t>& aln_iitree,
                      overlap_atomic_queue_t& ovlp_q,
                      range_atomic_queue_t& todo_in);

void write_graph_chunk(const seqindex_t& seqidx,
                       mmmulti::iitree<uint64_t, pos_t>& node_iitree,
                       mmmulti::iitree<uint64_t, pos_t>& path_iitree,
                       std::ofstream& seq_v_out,
                       std::map<pos_t, range_t>& range_buffer,
                       std::vector<std::pair<uint64_t, uint64_t>>* dsets_ptr,
                       uint64_t repeat_max,
                       uint64_t min_repeat_dist);

size_t compute_transitive_closures(
    const seqindex_t& seqidx,
    mmmulti::iitree<uint64_t, pos_t>& aln_iitree, // input alignment matches between query seqs
    const std::string& seq_v_file,
    mmmulti::iitree<uint64_t, pos_t>& node_iitree, // maps graph to input
    mmmulti::iitree<uint64_t, pos_t>& path_iitree, // maps input to graph
    uint64_t repeat_max,
    uint64_t min_repeat_dist,
    uint64_t transclose_batch_size,
    bool show_progress,
    uint64_t num_threads,
    const std::chrono::time_point<std::chrono::steady_clock>& start_time);

}
