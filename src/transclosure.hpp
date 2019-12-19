#ifndef TRANSCLOSURE_HPP_INCLUDED
#define TRANSCLOSURE_HPP_INCLUDED

#include <string>
#include <fstream>
#include <iostream>
#include <unordered_set>
#include <set>
#include "sdsl/bit_vectors.hpp"
#include "seqindex.hpp"
#include "mmiitree.hpp"
#include "pos.hpp"
#include "match.hpp"

namespace seqwish {


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
