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
#include "pos.hpp"

namespace seqwish {

using mmmultimap::multimap;

size_t compute_transitive_closures(
    seqindex_t& seqidx,
    multimap<uint64_t, pos_t>& aln_mm,
    const std::string& seq_v_file,
    multimap<uint64_t, uint64_t>& node_mm,
    multimap<uint64_t, uint64_t>& path_mm,
    uint64_t repeat_max);

}

#endif
