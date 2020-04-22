#include "compact.hpp"
#include "mmmultimap.hpp"

namespace seqwish {

void compact_nodes(
    seqindex_t& seqidx,
    size_t graph_size,
    mmmulti::iitree<uint64_t, pos_t>& node_iitree,
    mmmulti::iitree<uint64_t, pos_t>& path_iitree,
    sdsl::bit_vector& seq_id_bv) {
    // for each pair of positions in the graph base seq
    // do we have any links from the first that don't go to the second?
    // do we have any links to the second that don't come from the first?
    seq_id_bv[0] = 1; // set first node start
    size_t num_seqs = seqidx.n_seqs();
#pragma omp parallel for
    for (size_t i = 1; i <= num_seqs; ++i) {
        size_t j = seqidx.nth_seq_offset(i);
        size_t seq_len = seqidx.nth_seq_length(i);
        size_t k = j + seq_len;
        //std::cerr << "compact " << seqidx.nth_name(i) << " " << seqidx.nth_seq_length(i) << " " << j << " " << k << std::endl;
        while (j < k) {
            std::vector<size_t> ovlp;
            path_iitree.overlap(j, j+1, ovlp);
            // each input base should only map one place in the graph
            if (ovlp.size() != 1) {
                std::cerr << "found " << ovlp.size()  << " overlaps for seq " << seqidx.nth_name(i) << " idx " << i << " at j=" << j << " of " << k << std::endl;
                for (auto& o : ovlp) {
                    std::cerr << "ovlp_start_in_q = " << path_iitree.start(o) << " "
                              << "ovlp_end_in_q = " << path_iitree.end(o) << " "
                              << "pos_start_in_s = " << pos_to_string(path_iitree.data(o)) << std::endl;
                }
                assert(false);
            }
            size_t idx = ovlp.front();
            uint64_t ovlp_start_in_q = path_iitree.start(idx);
            uint64_t ovlp_end_in_q = path_iitree.end(idx);
            pos_t pos_start_in_s = path_iitree.data(idx);
            bool match_is_rev = is_rev(pos_start_in_s);
            // mark a node start and end
            pos_t pos_end_in_s = pos_start_in_s;
            if (!match_is_rev) {
                incr_pos(pos_end_in_s, ovlp_end_in_q - ovlp_start_in_q);
#pragma omp critical (seq_id_bv)
                {
                    //std::cerr << "marking node+ start " << offset(pos_start_in_s) << " of " << seq_id_bv.size() << std::endl;
                    seq_id_bv[offset(pos_start_in_s)] = 1;
                }
#pragma omp critical (seq_id_bv)
                {
                    //std::cerr << "marking node+ end " << offset(pos_end_in_s) << " of " << seq_id_bv.size() << std::endl;
                    seq_id_bv[offset(pos_end_in_s)] = 1;
                }
            } else {
                incr_pos(pos_end_in_s, ovlp_end_in_q - ovlp_start_in_q - 1);
#pragma omp critical (seq_id_bv)
                {
                    //std::cerr << "marking node- start " << offset(pos_end_in_s) << " of " << seq_id_bv.size() << std::endl;
                    seq_id_bv[offset(pos_end_in_s)] = 1;
                }
#pragma omp critical (seq_id_bv)
                {
                    //std::cerr << "marking node- end " << offset(pos_start_in_s) << " of " << seq_id_bv.size() << std::endl;
                    seq_id_bv[offset(pos_start_in_s)+1] = 1;
                }
            }
            j = ovlp_end_in_q;
        }
    }
    //std::cerr << graph_size << " " << seq_id_bv.size() << std::endl;
    seq_id_bv[graph_size] = 1;
}

}
