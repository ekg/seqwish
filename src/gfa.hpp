#ifndef GFA_HPP_INCLUDED
#define GFA_HPP_INCLUDED

#include <iostream>
#include <sstream>
#include "mmmultimap.hpp"
#include "seqindex.hpp"
#include "pos.hpp"
#include "mmap.hpp"

namespace seqwish {

using mmmultimap::multimap;

void emit_gfa(std::ostream& out,
              size_t graph_length,
              const std::string& seq_v_file,
              multimap<uint64_t, pos_t>& path_mm,
              multimap<pos_t, pos_t>& link_fwd_mm,
              multimap<pos_t, pos_t>& link_rev_mm,
              const sdsl::sd_vector<>& seq_id_cbv,
              const sdsl::sd_vector<>::rank_1_type& seq_id_cbv_rank,
              const sdsl::sd_vector<>::select_1_type& seq_id_cbv_select,
              seqindex_t& seqidx);

}

#endif
