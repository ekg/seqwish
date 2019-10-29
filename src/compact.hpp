#ifndef COMPACT_HPP_INCLUDED
#define COMPACT_HPP_INCLUDED

#include <vector>
#include "sdsl/bit_vectors.hpp"
#include "seqindex.hpp"
#include "mmiitree.hpp"
#include "iitii_types.hpp"
#include "pos.hpp"

namespace seqwish {


void compact_nodes(
    seqindex_t& seqidx,
    size_t graph_size,
    range_pos_iitii& node_iitree,
    range_pos_iitii& path_iitree,
    sdsl::bit_vector& seq_id_bv);

}

#endif
