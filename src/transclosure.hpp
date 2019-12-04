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
#include "iitii_types.hpp"
#include "pos.hpp"

namespace seqwish {


size_t compute_transitive_closures(
    seqindex_t& seqidx,
    range_pos_iitii& aln_iitree, // input alignment matches between query seqs
    const std::string& seq_v_file,
    mmmulti::iitree<uint64_t, pos_t>& node_iitree, // maps graph to input
    mmmulti::iitree<uint64_t, pos_t>& path_iitree, // maps input to graph
    uint64_t repeat_max,
    uint64_t min_transclose_len,
    uint64_t transclose_batch_size);

}

#endif
