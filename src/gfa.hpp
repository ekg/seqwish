#ifndef GFA_HPP_INCLUDED
#define GFA_HPP_INCLUDED

#include <iostream>
#include <sstream>
#include "mmiitree.hpp"
#include "seqindex.hpp"
#include "pos.hpp"
#include "mmap.hpp"

namespace seqwish {


void emit_gfa(std::ostream& out,
              size_t graph_length,
              const std::string& seq_v_file,
              mmmulti::iitree<uint64_t, pos_t>& node_iitree,
              mmmulti::iitree<uint64_t, pos_t>& path_iitree,
              const sdsl::sd_vector<>& seq_id_cbv,
              const sdsl::sd_vector<>::rank_1_type& seq_id_cbv_rank,
              const sdsl::sd_vector<>::select_1_type& seq_id_cbv_select,
              seqindex_t& seqidx,
              mmmulti::set<std::tuple<uint64_t, bool, uint64_t, bool>>& link_mmset);

}

#endif
