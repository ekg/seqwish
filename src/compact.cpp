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
            uint64_t overlap_count = 0;
            uint64_t ovlp_start_in_q;
            uint64_t ovlp_end_in_q;
            pos_t pos_start_in_s;
            path_iitree.overlap(
                j, j+1,
                [&](const uint64_t& start,
                    const uint64_t& end,
                    const pos_t& pos) {
                    ++overlap_count;
                    ovlp_start_in_q = start;
                    ovlp_end_in_q = end;
                    pos_start_in_s = pos;
                });
            // each input base should only map one place in the graph
            if (overlap_count != 1) {
                std::cerr << "[seqwish::compact] error: found " << overlap_count << " overlaps for seq " << seqidx.nth_name(i) << " idx " << i << " at j=" << j << " of " << k << std::endl;
                path_iitree.overlap(
                    j, j+1,
                    [&](const uint64_t& start,
                    const uint64_t& end,
                    const pos_t& pos) {
                        std::cerr << "ovlp_start_in_q = " << start << " "
                                  << "ovlp_end_in_q = " << end << " "
                                  << "pos_start_in_s = " << pos_to_string(pos) << std::endl;
                    });
                assert(false);
                exit(1);
            }

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
