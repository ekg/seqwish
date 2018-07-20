#ifndef TRANSCLOSURE_HPP_INCLUDED
#define TRANSCLOSURE_HPP_INCLUDED

#include <string>
#include <fstream>
#include <iostream>
#include <unordered_set>
#include <set>
#include "sdsl/bit_vectors.hpp"
#include "seqindex.hpp"
#include "dmultimap.hpp"
#include "pos.hpp"

namespace seqwish {

size_t compute_transitive_closures(
    seqindex_t& seqidx,
    dmultimap<uint64_t, aln_pos_t>& aln_mm,
    const std::string& seq_v_file,
    dmultimap<uint64_t, uint64_t>& node_mm,
    dmultimap<uint64_t, uint64_t>& path_mm,
    uint64_t repeat_max);

}

#endif
