#ifndef LINKS_HPP_INCLUDED
#define LINKS_HPP_INCLUDED

#include <vector>
#include <iostream>
#include "paryfor.hpp"
#include "seqindex.hpp"
#include "mmmultiset.hpp"
#include "mmiitree.hpp"
#include "pos.hpp"

namespace seqwish {


void derive_links(seqindex_t& seqidx,
                  mmmulti::iitree<uint64_t, pos_t>& node_iitree,
                  mmmulti::iitree<uint64_t, pos_t>& path_iitree,
                  const sdsl::sd_vector<>& seq_id_cbv,
                  const sdsl::sd_vector<>::rank_1_type& seq_id_cbv_rank,
                  const sdsl::sd_vector<>::select_1_type& seq_id_cbv_select,
                  mmmulti::set<std::pair<pos_t, pos_t>>& link_mmset,
                  const uint64_t& num_threads);

}

#endif
