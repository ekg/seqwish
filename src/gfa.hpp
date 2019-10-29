#ifndef GFA_HPP_INCLUDED
#define GFA_HPP_INCLUDED

#include <iostream>
#include <sstream>
#include "mmiitree.hpp"
#include "mmmultiset.hpp"
#include "seqindex.hpp"
#include "pos.hpp"
#include "mmap.hpp"
#include "iitii_types.hpp"

namespace seqwish {


void emit_gfa(std::ostream& out,
              size_t graph_length,
              const std::string& seq_v_file,
              range_pos_iitii& node_iitree,
              range_pos_iitii& path_iitree,
              const sdsl::sd_vector<>& seq_id_cbv,
              const sdsl::sd_vector<>::rank_1_type& seq_id_cbv_rank,
              const sdsl::sd_vector<>::select_1_type& seq_id_cbv_select,
              seqindex_t& seqidx,
              mmmulti::set<std::pair<pos_t, pos_t>>& link_mmset);

}

#endif
