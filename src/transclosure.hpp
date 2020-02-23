#pragma once

#include <string>
#include <fstream>
#include <iostream>
#include <unordered_set>
#include <set>
#include <thread>
#include "sdsl/bit_vectors.hpp"
#include "atomic_bitvector.hpp"
#include "seqindex.hpp"
#include "mmiitree.hpp"
#include "pos.hpp"
#include "match.hpp"
#include "ips4o.hpp"
#include "spinlock.hpp"
#include "dset64-gccAtomic.hpp"
#include "atomic_queue.h"

namespace seqwish {

typedef atomic_queue::AtomicQueue2<std::pair<pos_t, uint64_t>, 2 << 16> range_atomic_queue_t;

struct range_t {
    //uint64_t seq_id = 0;
    uint64_t begin = 0;
    uint64_t end = 0;
};

void extend_range(const uint64_t& s_pos,
                  const pos_t& q_pos,
                  std::map<pos_t, range_t>& range_buffer);

void flush_ranges(const uint64_t& s_pos,
                  std::map<pos_t, range_t>& range_buffer,
                  mmmulti::iitree<uint64_t, pos_t>& node_iitree,
                  mmmulti::iitree<uint64_t, pos_t>& path_iitree);

void for_each_fresh_range(const match_t& range,
                          atomicbitvector::atomic_bv_t& seen_bv,
                          const seqindex_t& seqidx,
                          const std::function<void(match_t)>& lambda);

void handle_range(match_t s,
                  atomicbitvector::atomic_bv_t& seen_bv,
                  atomicbitvector::atomic_bv_t& curr_bv,
                  const seqindex_t& seqidx,
                  const uint64_t& query_start,
                  const uint64_t& query_end,
                  std::vector<std::pair<match_t, bool>>& ovlp,
                  range_atomic_queue_t& todo,
                  std::vector<std::pair<pos_t, uint64_t>>& overflow);

void explore_overlaps(const match_t& b,
                      atomicbitvector::atomic_bv_t& seen_bv,
                      atomicbitvector::atomic_bv_t& curr_bv,
                      const seqindex_t& seqidx,
                      mmmulti::iitree<uint64_t, pos_t>& aln_iitree,
                      std::vector<std::pair<match_t, bool>>& ovlp,
                      range_atomic_queue_t& todo,
                      std::vector<std::pair<pos_t, uint64_t>>& overflow);

size_t compute_transitive_closures(
    const seqindex_t& seqidx,
    mmmulti::iitree<uint64_t, pos_t>& aln_iitree, // input alignment matches between query seqs
    const std::string& seq_v_file,
    mmmulti::iitree<uint64_t, pos_t>& node_iitree, // maps graph to input
    mmmulti::iitree<uint64_t, pos_t>& path_iitree, // maps input to graph
    uint64_t repeat_max,
    uint64_t transclose_batch_size);

}
