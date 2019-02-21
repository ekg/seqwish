#include "vgp.hpp"
#include "dmultimap.hpp"

namespace seqwish {

void emit_vgp(const std::string& basename,
              size_t graph_length,
              const std::string& seq_v_file,
              dmultimap<uint64_t, pos_t>& path_mm,
              dmultimap<pos_t, pos_t>& link_fwd_mm,
              dmultimap<pos_t, pos_t>& link_rev_mm,
              const sdsl::sd_vector<>& seq_id_cbv,
              const sdsl::sd_vector<>::rank_1_type& seq_id_cbv_rank,
              const sdsl::sd_vector<>::select_1_type& seq_id_cbv_select,
              seqindex_t& seqidx) {

    // the contig or scaffold sequences, which are the nodes in our graph
    std::string seq_file = basename+".seq";
    std::ofstream seq_out(seq_file);

    // the joins described by the graph
    std::string scf_file = basename+".scf";
    std::ofstream scf_out(scf_file);

    // the mapping between the input sequences and the output nodes/contigs/scaffolds
    std::string sxs_file = basename+".sxs";
    std::ofstream sxs_out(sxs_file);    
    
    //out << "H" << "\t" << "VN:Z:1.0" << std::endl;
    int seq_v_fd = -1;
    char* seq_v_buf = nullptr;
    size_t seq_v_filesize = mmap_open(seq_v_file, seq_v_buf, seq_v_fd);
    // write the nodes
    // these are delimited in the seq_v_file by the markers in seq_id_civ
    auto show_links = [&](const pos_t& p) { std::cerr << pos_to_string(p) << " " << pos_to_string(make_pos_t(seq_id_cbv_rank(offset(p)), is_rev(p))) << ", "; };
    size_t n_nodes = seq_id_cbv_rank(seq_id_cbv.size());
    for (size_t id = 1; id <= n_nodes; ++id) {
        size_t node_start = seq_id_cbv_select(id);
        size_t node_length = (id==n_nodes ? seq_id_cbv.size() : seq_id_cbv_select(id+1)) - node_start;
        //std::cerr << id << " "  << node_start << " " << node_length << std::endl;
        char seqc[node_length];
        memcpy(&seqc[0], &seq_v_buf[node_start], node_length);
        std::string seq(seqc, node_length);
        seq_out << "S" << "\t" << seq.size() << "\t" << seq << std::endl;
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

        auto print_to_join = [&scf_out, &id, &seq_id_cbv, &seq_id_cbv_rank](const pos_t& p) {
            size_t i = offset(p);
            // internal links which do not go to the head or tail of a compacted node
            if (!is_rev(p) && (seq_id_cbv[i]||i==seq_id_cbv.size())
                || is_rev(p) && seq_id_cbv[i-1]) {
                pos_t node = make_pos_t(seq_id_cbv_rank(i), is_rev(p));
                scf_out << "J" << "\t"
                        << offset(node) << "\t"
                        << id << "\t"
                        << (is_rev(node)?"s":"e") << "\t"
                        << "s" // forward relative
                        << std::endl;
            }
        };

        auto print_from_join = [&scf_out, &id, &seq_id_cbv, &seq_id_cbv_rank](const pos_t& p) {
            size_t i = offset(p);
            // internal links which do not go to the head or tail of a compacted node
            if (!is_rev(p) && seq_id_cbv[i-1]
                || is_rev(p) && (seq_id_cbv[i] || i == seq_id_cbv.size())) {
                pos_t node = make_pos_t(seq_id_cbv_rank(i), is_rev(p));
                scf_out << "J" << "\t"
                        << id << "\t" 
                        << offset(node) << "\t"
                        << "e" << "\t"
                        << (is_rev(node)?"e":"s")
                        << std::endl;
            }
        };
        
        //std::cerr << "fwd start" << std::endl;
        link_rev_mm.for_unique_values_of(node_start_fwd, print_to_join);
        //std::cerr << "fwd end" << std::endl;
        link_fwd_mm.for_unique_values_of(node_end_fwd, print_from_join);

    }

    // write the paths
    // iterate over the sequence positions, emitting a node at every edge crossing
    size_t num_seqs = seqidx.n_seqs();
    for (size_t i = 1; i <= num_seqs; ++i) {
        size_t j = seqidx.nth_seq_offset(i) + 1;
        size_t k = j+seqidx.nth_seq_length(i);
        //std::cerr << seqidx.nth_name(i) << " " << seqidx.nth_seq_length(i) << " " << j << " " << k << std::endl;
        std::vector<pos_t> path_v;
        pos_t last_pos = 0;
        pos_t last_node = 0;
        for ( ; j < k; ++j) {
            std::vector<pos_t> v = path_mm.values(j);
            // each input base should only map one place in the graph
            assert(v.size() == 1);
            auto& p = v.front();
            //std::cerr << j << " -> " << pos_to_string(p) << std::endl;
            // validate the path
            char c = seq_v_buf[offset(p)-1];
            if (is_rev(p)) c = dna_reverse_complement(c);
            //std::cerr << seqidx.at_pos(make_pos_t(j, false)) << " " << c << std::endl;
            assert(seqidx.at_pos(make_pos_t(j, false)) == c);
            pos_t node = make_pos_t(seq_id_cbv_rank(offset(p)), is_rev(p));
            pos_t lp = last_pos; incr_pos(lp);
            //std::cerr << offset(last_node) << std::endl;
            if (offset(last_node)
                && (p != lp
                    || node != last_node)) { // or if we
                //out << pos_to_string(last_node) << ",";
                path_v.push_back(last_node);
            }
            last_pos = p;
            last_node = node;
        }
        path_v.push_back(last_node);
        std::stringstream cigarss;
        std::stringstream pathss;
        uint64_t pos_in_b = 0;
        for (auto& p : path_v) {
            uint64_t id = offset(p);
            size_t node_start = seq_id_cbv_select(id);
            size_t node_length = (id==n_nodes ? seq_id_cbv.size() : seq_id_cbv_select(id+1)) - node_start;
            uint64_t start_pos = (is_rev(p) ? pos_in_b+node_length : pos_in_b);
            uint64_t end_pos = (is_rev(p) ? pos_in_b : pos_in_b+node_length);
            //cigarss << node_length << "M" << ",";
            sxs_out << "A" << "\t"
                    << id << "\t"
                    << seqidx.nth_name(i) << "\n"
                    << "I" << "\t"
                    << 0 << "\t" << node_length << "\t"
                    << start_pos << "\t" << end_pos << "\n";
            
            pos_in_b += node_length;
        }
        //pathss.seekp(-1, pathss.cur); pathss << '\t';
        //cigarss.seekp(-1, cigarss.cur); cigarss << std::endl;
    }

    mmap_close(seq_v_buf, seq_v_fd, seq_v_filesize);

}

}
