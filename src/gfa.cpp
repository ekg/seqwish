#include "gfa.hpp"
#include "dmultimap.hpp"

namespace seqwish {

void emit_gfa(std::ostream& out,
              size_t graph_length,
              const std::string& seq_v_file,
              dmultimap<uint64_t, pos_t>& path_mm,
              dmultimap<pos_t, pos_t>& link_fwd_mm,
              dmultimap<pos_t, pos_t>& link_rev_mm,
              const sdsl::sd_vector<>& seq_id_cbv,
              const sdsl::sd_vector<>::rank_1_type& seq_id_cbv_rank,
              const sdsl::sd_vector<>::select_1_type& seq_id_cbv_select) {

    out << "H" << "\t" << "VN:Z:1.0" << std::endl;
    std::ifstream seq_in(seq_v_file.c_str());
    // write the nodes
    // these are delimited in the seq_v_file by the markers in seq_id_civ
    auto show_links = [&](const pos_t& p) { std::cerr << pos_to_string(p) << " " << pos_to_string(make_pos_t(seq_id_cbv_rank(offset(p)), is_rev(p))) << ", "; };
    size_t n_nodes = seq_id_cbv_rank(seq_id_cbv.size());
    for (size_t id = 1; id <= n_nodes; ++id) {
        size_t node_start = seq_id_cbv_select(id);
        size_t node_length = (id==n_nodes ? seq_id_cbv.size() : seq_id_cbv_select(id+1)) - node_start;
        //std::cerr << id << " "  << node_start << " " << node_length << std::endl;
        char seq[node_length+1]; seq[node_length]='\0';
        seq_in.seekg(node_start);
        seq_in.read(seq, node_length);
        out << "S" << "\t" << id << "\t" << std::string(seq) << std::endl;
        // get the links of this node
        // to the forward or reverse start
        pos_t node_start_fwd = make_pos_t(node_start+1, false);
        pos_t node_end_fwd = make_pos_t(node_start+node_length, false);
        // from the forward or reverse end
        pos_t node_start_rev = make_pos_t(node_start+node_length, true);
        pos_t node_end_rev = make_pos_t(node_start+1, true);

        /*
        std::cerr << "node extents "
                  << pos_to_string(node_start_fwd) << " "
                  << pos_to_string(node_end_fwd) << " "
                  << pos_to_string(node_start_rev) << " "
                  << pos_to_string(node_end_rev) << std::endl;
        std::cerr << "things to this node fwd" << " ";
        link_rev_mm.for_values_of(node_start_fwd, show_links);
        std::cerr << std::endl;
        std::cerr << "things to this node rev" << " ";
        link_rev_mm.for_values_of(node_start_rev, show_links);
        std::cerr << std::endl;
        std::cerr << "things from this node fwd" << " ";
        link_fwd_mm.for_values_of(node_end_fwd, show_links);
        std::cerr << std::endl;
        std::cerr << "things from this node rev" << " ";
        link_fwd_mm.for_values_of(node_end_rev, show_links);
        std::cerr << std::endl;
        */

        auto print_to_link = [&out, &id, &seq_id_cbv, &seq_id_cbv_rank](const pos_t& p) {
            size_t i = offset(p);
            // internal links which do not go to the head or tail of a compacted node
            if (!is_rev(p) && (seq_id_cbv[i]||i==seq_id_cbv.size())
                || is_rev(p) && seq_id_cbv[i-1]) {
                pos_t node = make_pos_t(seq_id_cbv_rank(i), is_rev(p));
                out << "L" << "\t"
                    << offset(node) << "\t" << (is_rev(node)?"-":"+") << "\t"
                    << id << "\t" << "+" << "\t"
                    << "OM" << std::endl;
            }
        };

        auto print_from_link = [&out, &id, &seq_id_cbv, &seq_id_cbv_rank](const pos_t& p) {
            size_t i = offset(p);
            // internal links which do not go to the head or tail of a compacted node
            if (!is_rev(p) && seq_id_cbv[i-1]
                || is_rev(p) && (seq_id_cbv[i] || i == seq_id_cbv.size())) {
                pos_t node = make_pos_t(seq_id_cbv_rank(i), is_rev(p));
                out << "L" << "\t"
                    << id << "\t" << "+" << "\t"
                    << offset(node) << "\t" << (is_rev(node)?"-":"+") << "\t"
                    << "OM" << std::endl;
            }
        };
        
        //std::cerr << "fwd start" << std::endl;
        link_rev_mm.for_unique_values_of(node_start_fwd, print_to_link);
        //std::cerr << "fwd end" << std::endl;
        link_fwd_mm.for_unique_values_of(node_end_fwd, print_from_link);
        
    }
    // write the paths
    
}

}
