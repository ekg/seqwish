#include "gfa.hpp"

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
              mmmulti::set<std::pair<pos_t, pos_t>>& link_mmset) {

    out << "H" << "\t" << "VN:Z:1.0" << std::endl;
    int seq_v_fd = -1;
    char* seq_v_buf = nullptr;
    size_t seq_v_filesize = mmap_open(seq_v_file, seq_v_buf, seq_v_fd);
    // write the nodes
    // these are delimited in the seq_v_file by the markers in seq_id_civ
    auto show_links = [&](const pos_t& p) { std::cerr << pos_to_string(p) << " " << pos_to_string(make_pos_t(seq_id_cbv_rank(offset(p)), is_rev(p))) << ", "; };
    size_t n_nodes = seq_id_cbv_rank(seq_id_cbv.size()-1);
    for (size_t id = 1; id <= n_nodes; ++id) {
        //std::cerr << "id " << id << " n_nodes " << n_nodes << std::endl;
        size_t node_start = seq_id_cbv_select(id);
        //size_t node_length = (id==n_nodes ? seq_id_cbv.size() : seq_id_cbv_select(id+1)) - node_start;
        size_t node_length = seq_id_cbv_select(id+1) - node_start;
        //std::cerr << id << " "  << node_start << " " << node_length << std::endl;
        std::string seq; seq.resize(node_length);
        memcpy((void*)seq.c_str(), &seq_v_buf[node_start], node_length);
        out << "S" << "\t" << id << "\t" << seq << std::endl;
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

    }

    auto print_link = [&out](const std::pair<pos_t, pos_t>& p) {
        auto& from = p.first;
        auto& to = p.second;
        if (from && to) {
            out << "L" << "\t"
                << offset(from) << "\t" << (is_rev(from)?"-":"+") << "\t"
                << offset(to) << "\t" << (is_rev(to)?"-":"+") << "\t"
                << "0M" << std::endl;
        }
    };
    link_mmset.for_each_unique_value(print_link);

    // write the paths
    // iterate over the sequence positions, emitting a node at every edge crossing
    size_t num_seqs = seqidx.n_seqs();
    for (size_t i = 1; i <= num_seqs; ++i) {
        size_t j = seqidx.nth_seq_offset(i)+1;
        size_t seq_len = seqidx.nth_seq_length(i);
        size_t k = j + seq_len;
        //std::cerr << seqidx.nth_name(i) << " " << seqidx.nth_seq_length(i) << " " << j << " " << k << std::endl;
        std::vector<pos_t> path_v;
        uint64_t seen_bp = 0;
        uint64_t accumulated_bp = 0;
        while (j < k) {
            std::vector<size_t> ovlp;
            path_iitree.overlap(j, j+1, ovlp);
            // each input base should only map one place in the graph
            assert(ovlp.size() == 1);
            size_t idx = ovlp.front();
            uint64_t ovlp_start_in_q = path_iitree.start(idx);
            uint64_t ovlp_end_in_q = path_iitree.end(idx);
            pos_t pos_start_in_s = path_iitree.data(idx);
            bool match_is_rev = is_rev(pos_start_in_s);
            // iterate through the nodes in this range
            uint64_t length = ovlp_end_in_q - ovlp_start_in_q;
            //std::cerr << "overlap in path_iitree " << ovlp_start_in_q << ".." << ovlp_end_in_q << " " << pos_to_string(pos_start_in_s) << std::endl;
            pos_t q = make_pos_t(ovlp_start_in_q, false);
            pos_t p = pos_start_in_s;
            // validate the path
            // for each base in the range pointed to by the match, check that the input sequence we're processing matches the graph
            for (uint64_t k = 0; k < length; ++k) {
                if (seq_id_cbv[offset(p)-1]) {
                    uint64_t node_id = seq_id_cbv_rank(offset(p));
                    //std::cerr << "got to node " << node_id << (match_is_rev ? "-" : "+") << std::endl;
                    path_v.push_back(make_pos_t(node_id, match_is_rev));
                }
                /*else if (match_is_rev && seq_id_cbv[offset(p)]) {
                    uint64_t node_id = seq_id_cbv_rank(offset(p)-1);
                    std::cerr << "got to node " << node_id << "-" << std::endl;
                    path_v.push_back(make_pos_t(node_id, match_is_rev));
                }
                */
                char c = seq_v_buf[offset(p)-1];
                if (is_rev(p)) c = dna_reverse_complement(c);
                //std::cerr << pos_to_string(q) << " -> " << pos_to_string(p) << " " << seqidx.at_pos(q) << " vs " << c << std::endl;
                if (seqidx.at_pos(q) != c) {
                    std::cerr << "GRAPH BROKEN @ "
                        << seqidx.nth_name(i) << " " << pos_to_string(q) << " -> "
                        << pos_to_string(q) << std::endl;
                    assert(false);
                    exit(1); // for release builds
                }
                incr_pos(p, 1);
                incr_pos(q, 1);
            }
            seen_bp += length;
            j = ovlp_end_in_q;
        }
        if (seen_bp != seq_len) {
            std::cerr << "length for " << seqidx.nth_name(i) << ", expected " << seqidx.nth_seq_length(i) << " but got " << accumulated_bp << std::endl;
            assert(false);
            exit(1); // for release builds
        }
        std::stringstream cigarss;
        std::stringstream pathss;

        for (auto& p : path_v) {
            pathss << pos_to_string(p) << ",";
        }
        pathss.seekp(-1, pathss.cur); // trim the last ","
        cigarss << "*";
        if (path_v.size() > 2) {
            for (uint64_t q = 0; q < path_v.size()-2; ++q) {
                cigarss << ",*";
            }
        }
        pathss << '\t';
        //cigarss.seekp(-1, cigarss.cur);
        cigarss << std::endl;
        out << "P" << "\t"
            << seqidx.nth_name(i) << "\t"
            << pathss.str()
            << cigarss.str();
    }

    mmap_close(seq_v_buf, seq_v_fd, seq_v_filesize);

}

}
