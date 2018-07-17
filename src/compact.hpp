#ifndef COMPACT_HPP_INCLUDED
#define COMPACT_HPP_INCLUDED

#include <vector>
#include "sdsl/bit_vectors.hpp"
#include "seqindex.hpp"
#include "dmultimap.hpp"
#include "pos.hpp"

namespace seqwish {

void compact_nodes(
    size_t graph_size,
    dmultimap<pos_t, pos_t>& link_fwd_mm,
    dmultimap<pos_t, pos_t>& link_rev_mm,
    sdsl::bit_vector& seq_id_bv);

}

#endif
