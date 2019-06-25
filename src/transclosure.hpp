#ifndef TRANSCLOSURE_HPP_INCLUDED
#define TRANSCLOSURE_HPP_INCLUDED

#include <string>
#include <fstream>
#include <iostream>
#include <unordered_set>
#include <set>
#include "sdsl/bit_vectors.hpp"
#include "seqindex.hpp"
#include "mmmultimap.hpp"
#include "mmiitree.hpp"
#include "pos.hpp"

namespace seqwish {


size_t compute_transitive_closures(
    seqindex_t& seqidx,
    mmmulti::iitree<uint64_t, pos_t>& aln_iitree,
    const std::string& seq_v_file,
    mmmulti::map<uint64_t, uint64_t>& node_mm,
    mmmulti::map<uint64_t, uint64_t>& path_mm,
    uint64_t repeat_max,
    uint64_t min_transclose_len);

}

#endif
