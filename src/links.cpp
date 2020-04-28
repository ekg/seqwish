#include "links.hpp"

namespace seqwish {

void derive_links(seqindex_t& seqidx,
                  mmmulti::iitree<uint64_t, pos_t>& node_iitree,
                  mmmulti::iitree<uint64_t, pos_t>& path_iitree,
                  const sdsl::sd_vector<>& seq_id_cbv,
                  const sdsl::sd_vector<>::rank_1_type& seq_id_cbv_rank,
                  const sdsl::sd_vector<>::select_1_type& seq_id_cbv_select,
                  mmmulti::set<std::pair<pos_t, pos_t>>& link_mmset,
                  const uint64_t& num_threads) {
    // for each marked node
    // determine our edge context using the node_iitree and path_iitree
    // and write it into the mmset
    link_mmset.open_writer();
    size_t n_nodes = seq_id_cbv_rank(seq_id_cbv.size()-1);
//#pragma omp parallel for
    for (size_t id = 1; id <= n_nodes; ++id) {
        uint64_t node_start_in_s = seq_id_cbv_select(id); // select is 1-based
        uint64_t node_end_in_s = seq_id_cbv_select(id+1);
        //std::cerr << "links for node " << id << " start " << node_start_in_s << " end " << node_end_in_s << std::endl;
        // find the things on both sides of our node by looking in the node_iitree, finding what bits of the paths (in Q)
        // are there, and seeing what's on either side of them to decide what links we need
        node_iitree.overlap(
            node_start_in_s, node_end_in_s,
            [&](const uint64_t& ovlp_start_in_s,
                const uint64_t& ovlp_end_in_s,
                const pos_t& base_pos_start_in_q) {
                uint64_t ovlp_length = ovlp_end_in_s - ovlp_start_in_s;
                pos_t pos_start_in_q = base_pos_start_in_q;
                pos_t pos_end_in_q = pos_start_in_q;
                incr_pos(pos_end_in_q, ovlp_length - (ovlp_end_in_s - node_end_in_s));
                incr_pos(pos_start_in_q, node_start_in_s - ovlp_start_in_s);
                // get the positions in s on either side of this range in Q
                // note the need to handle orientation
                // determine which sequence we're in
                // find its boundaries
                bool curr_step_is_rev = is_rev(pos_start_in_q);
                uint64_t start_in_q = (curr_step_is_rev ? offset(pos_end_in_q)+1 : offset(pos_start_in_q));
                uint64_t end_in_q = (curr_step_is_rev ? offset(pos_start_in_q)+1 : offset(pos_end_in_q));
                uint64_t seq_id = seqidx.seq_id_at(start_in_q);
                uint64_t id_end = seqidx.seq_id_at(end_in_q-1);
                assert(seq_id == id_end);
                uint64_t seq_start = seqidx.nth_seq_offset(seq_id);
                uint64_t seq_end = seq_start + seqidx.nth_seq_length(seq_id);
                // and only consider cases where we'd be within the boundaries
                if (end_in_q+1 <= seq_end) {
                    path_iitree.overlap(
                        end_in_q, end_in_q+1,
                        [&](const uint64_t& ovlp_start_in_q,
                            const uint64_t& ovlp_end_in_q,
                            const pos_t& base_pos_start_in_s) {
                            pos_t pos_start_in_s = base_pos_start_in_s;
                            if (ovlp_start_in_q < end_in_q) {
                                //std::cerr << "contained in previous range on q" << std::endl;
                                incr_pos(pos_start_in_s, end_in_q - ovlp_start_in_q);
                            }
                            // find which node we're in here, and record a link
                            uint64_t next_id = seq_id_cbv_rank(offset(pos_start_in_s)+1);
                            bool next_step_is_rev = is_rev(pos_start_in_s);
                            link_mmset.append(std::make_pair(make_pos_t(id, curr_step_is_rev), make_pos_t(next_id, next_step_is_rev)));
                        });
                }
            });
    }
    link_mmset.index(num_threads);
}

}
