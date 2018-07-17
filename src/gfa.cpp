#include "gfa.hpp"
#include "dmultimap.hpp"

namespace seqwish {

void emit_gfa(std::ostream& out,
              size_t graph_length,
              const std::string& seq_v_file,
              dmultimap<uint64_t, pos_t>& path_mm,
              dmultimap<pos_t, pos_t>& link_fwd_mm,
              dmultimap<pos_t, pos_t>& link_rev_mm,
              const sdsl::sd_vector<>& seq_id_civ,
              const sdsl::sd_vector<>::rank_1_type& seq_id_cbv_rank,
              const sdsl::sd_vector<>::select_1_type& seq_id_cbv_select) {

    out << "H" << "\t" << "VN:Z:1.0" << std::endl;
    // write the nodes
    // write the links
    // write the paths
    
}

}
