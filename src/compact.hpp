#pragma once

#include <vector>
#include "sdsl/bit_vectors.hpp"
#include "atomic_bitvector.hpp"
#include "seqindex.hpp"
#include "mmiitree.hpp"
#include "pos.hpp"
#include "paryfor.hpp"

namespace seqwish {


void compact_nodes(
    seqindex_t& seqidx,
    size_t graph_size,
    mmmulti::iitree<uint64_t, pos_t>& node_iitree,
    mmmulti::iitree<uint64_t, pos_t>& path_iitree,
    sdsl::bit_vector& seq_id_bv,
    uint64_t num_threads);

}
