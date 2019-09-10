#include "gfa.hpp"
#include "mmmultimap.hpp"

namespace seqwish {

void emit_gfa(std::ostream& out,
              size_t graph_length,
              const std::string& seq_v_file,
              mmmulti::map<uint64_t, pos_t>& path_mm,
              mmmulti::map<pos_t, pos_t>& link_fwd_mm,
              mmmulti::map<pos_t, pos_t>& link_rev_mm,
              const sdsl::sd_vector<>& seq_id_cbv,
              const sdsl::sd_vector<>::rank_1_type& seq_id_cbv_rank,
              const sdsl::sd_vector<>::select_1_type& seq_id_cbv_select,
              seqindex_t& seqidx) {

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

        auto print_to_link = [&out, &id, &seq_id_cbv, &seq_id_cbv_rank](const pos_t& p) {
            size_t i = offset(p);
            // internal links which do not go to the head or tail of a compacted node
            if (!is_rev(p) && seq_id_cbv[i]
                || is_rev(p) && seq_id_cbv[i-1]) {
                pos_t node = make_pos_t(seq_id_cbv_rank(i), is_rev(p));
                out << "L" << "\t"
                    << offset(node) << "\t" << (is_rev(node)?"-":"+") << "\t"
                    << id << "\t" << "+" << "\t"
                    << "0M" << std::endl;
            }
        };

        auto print_from_link = [&out, &id, &seq_id_cbv, &seq_id_cbv_rank](const pos_t& p) {
            size_t i = offset(p);
            // internal links which do not go to the head or tail of a compacted node
            if (!is_rev(p) && seq_id_cbv[i-1]
                || is_rev(p) && seq_id_cbv[i]) {
                pos_t node = make_pos_t(seq_id_cbv_rank(i), is_rev(p));
                out << "L" << "\t"
                    << id << "\t" << "+" << "\t"
                    << offset(node) << "\t" << (is_rev(node)?"-":"+") << "\t"
                    << "0M" << std::endl;
            }
        };
        
        //std::cerr << "fwd start" << std::endl;
        link_rev_mm.for_unique_values_of(node_start_fwd, print_to_link);
        //std::cerr << "fwd end" << std::endl;
        link_fwd_mm.for_unique_values_of(node_end_fwd, print_from_link);
        //pos_t node_start_rev = make_pos_t(node_start+node_length, true);
        //pos_t node_end_rev = make_pos_t(node_start+1, true);

    }

    // write the paths
    // iterate over the sequence positions, emitting a node at every edge crossing
    size_t num_seqs = seqidx.n_seqs();
    for (size_t i = 1; i <= num_seqs; ++i) {
        size_t j = seqidx.nth_seq_offset(i);
        size_t seq_len = seqidx.nth_seq_length(i);
        size_t k = j + seq_len;
        //std::cerr << seqidx.nth_name(i) << " " << seqidx.nth_seq_length(i) << " " << j << " " << k << std::endl;
        std::vector<pos_t> path_v;
        uint64_t seen_bp = 0;
        uint64_t accumulated_bp = 0;
        for ( ; j < k; ++j) {
            std::vector<pos_t> v = path_mm.values(j+1);
            // each input base should only map one place in the graph
            assert(v.size() == 1);
            auto& p = v.front();
            // validate the path
            char c = seq_v_buf[offset(p)-1];
            if (is_rev(p)) c = dna_reverse_complement(c);
            assert(seqidx.at_pos(make_pos_t(j+1, false)) == c);
            // are we at the start of a node?
            // if so, write to the path
            //std::cerr << seqidx.nth_name(i) << " " << seqidx.nth_seq_length(i) << " @" << offset(p)-1 << " in seq_v, seen_bp " << seen_bp << " " << accumulated_bp << " " << c << " == " << seqidx.at_pos(make_pos_t(j+1, false)) << std::endl;
            if (!is_rev(p)) {
                if (seq_id_cbv[offset(p)-1] == 1) {
                    size_t id = seq_id_cbv_rank(offset(p));
                    path_v.push_back(make_pos_t(id, is_rev(p)));
                    //std::cerr << "adding " << id << "+" << std::endl;
                    // TODO why does this break in some cases?
                    //assert(seen_bp == accumulated_bp);
                    accumulated_bp += seq_id_cbv_select(id+1) - seq_id_cbv_select(id);
                }
            } else {
                if (seq_id_cbv[offset(p)] == 1) {
                    size_t id = seq_id_cbv_rank(offset(p));
                    path_v.push_back(make_pos_t(id, is_rev(p)));
                    //std::cerr << "adding " << id << "-" << std::endl;
                    // TODO why does this break in some cases?
                    //assert(seen_bp == accumulated_bp);
                    accumulated_bp += seq_id_cbv_select(id+1) - seq_id_cbv_select(id);
                }
            }
            ++seen_bp;
        }
        if (accumulated_bp != seq_len) {
            std::cerr << "length for " << seqidx.nth_name(i) << ", expected " << seqidx.nth_seq_length(i) << " but got " << accumulated_bp << std::endl;
            assert(false);
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
