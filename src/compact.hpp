#ifndef COMPACT_HPP_INCLUDED
#define COMPACT_HPP_INCLUDED

#include <vector>
#include "sdsl/bit_vectors.hpp"
#include "seqindex.hpp"
#include "mmmultimap.hpp"
#include "pos.hpp"

namespace seqwish {


void compact_nodes(
    seqindex_t& seqidx,
    size_t graph_size,
    mmmulti::map<uint64_t, pos_t>& path_mm,
    mmmulti::map<pos_t, pos_t>& link_fwd_mm,
    mmmulti::map<pos_t, pos_t>& link_rev_mm,
    sdsl::bit_vector& seq_id_bv);

}

#endif
