#ifndef TRANSCLOSURE_HPP_INCLUDED
#define TRANSCLOSURE_HPP_INCLUDED

#include <string>
#include <fstream>
#include <iostream>
#include <unordered_set>
#include <set>
#include "sdsl/bit_vectors.hpp"
#include "atomic_bitvector.hpp"
#include "seqindex.hpp"
#include "mmiitree.hpp"
#include "pos.hpp"
#include "match.hpp"
#include "ips4o.hpp"
#include "spinlock.hpp"
#include "dset64-gccAtomic.hpp"
#include "BooPHF.h"

namespace seqwish {

void extend_range(const uint64_t& s_pos,
                  const pos_t& q_pos,
                  std::map<pos_t, std::pair<uint64_t, uint64_t>>& range_buffer);

void flush_ranges(const uint64_t& s_pos,
                  std::map<pos_t, std::pair<uint64_t, uint64_t>>& range_buffer,
                  mmmulti::iitree<uint64_t, pos_t>& node_iitree,
                  mmmulti::iitree<uint64_t, pos_t>& path_iitree);

void for_each_fresh_range(const match_t& range,
                          atomicbitvector::atomic_bv_t& seen_bv,
                          const std::function<void(match_t)>& lambda);

void handle_range(match_t s,
                  atomicbitvector::atomic_bv_t& curr_bv,
                  const uint64_t& query_start,
                  const uint64_t& query_end,
                  std::vector<std::pair<match_t, bool>>& ovlp,
                  std::set<std::pair<pos_t, uint64_t>>& todo);

void explore_overlaps(const match_t& b,
                      atomicbitvector::atomic_bv_t& seen_bv,
                      atomicbitvector::atomic_bv_t& curr_bv,
                      mmmulti::iitree<uint64_t, pos_t>& aln_iitree,
                      std::vector<std::pair<match_t, bool>>& ovlp,
                      std::set<std::pair<pos_t, uint64_t>>& todo);

size_t compute_transitive_closures(
    seqindex_t& seqidx,
    mmmulti::iitree<uint64_t, pos_t>& aln_iitree, // input alignment matches between query seqs
    const std::string& seq_v_file,
    mmmulti::iitree<uint64_t, pos_t>& node_iitree, // maps graph to input
    mmmulti::iitree<uint64_t, pos_t>& path_iitree, // maps input to graph
    uint64_t repeat_max,
    uint64_t min_transclose_len,
    uint64_t transclose_batch_size);

class pos_t_hasher {
public:
    // the class should have operator () with this signature         
    uint64_t operator () (const uint64_t& key, uint64_t seed=0) const {
        uint64_t hash = hash_fn(key);
        hash ^= seed;
        return hash;
    }
    std::hash<uint64_t> hash_fn;
};

}

#endif
