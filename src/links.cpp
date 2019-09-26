#include "links.hpp"

namespace seqwish {

void derive_links(seqindex_t& seqidx,
                  mmmulti::iitree<uint64_t, pos_t>& node_iitree,
                  mmmulti::iitree<uint64_t, pos_t>& path_iitree,
                  const sdsl::sd_vector<>& seq_id_cbv,
                  const sdsl::sd_vector<>::rank_1_type& seq_id_cbv_rank,
                  const sdsl::sd_vector<>::select_1_type& seq_id_cbv_select,
                  mmmulti::set<std::pair<pos_t, pos_t>>& link_mmset) {
    // for each marked node
    // determine our edge context using the node_iitree and path_iitree
    // and write it into the mmset
    size_t n_nodes = seq_id_cbv_rank(seq_id_cbv.size()-1);
    for (size_t id = 1; id <= n_nodes; ++id) {
        size_t node_start = seq_id_cbv_select(id)+1; // node_iitree and path_iitree are 1-based
        size_t node_end = seq_id_cbv_select(id+1)+1;
        std::cerr << "node " << id << " start " << node_start << " end " << node_end << std::endl;
        // find the things on both sides of our node by looking in the node_iitree
        std::vector<size_t> node_ovlp;
        node_iitree.overlap(node_start, node_end, node_ovlp);
        for (auto& idx : node_ovlp) {
            uint64_t start_in_s = node_iitree.start(idx);
            uint64_t end_in_s = node_iitree.end(idx);
            uint64_t length = end_in_s - start_in_s;
            pos_t pos_start_in_q = node_iitree.data(idx);
            pos_t pos_end_in_q = pos_start_in_q;
            incr_pos(pos_end_in_q, length - (end_in_s - node_end));
            incr_pos(pos_start_in_q, node_start - start_in_s);
            length -= (node_start - start_in_s) + (node_end - end_in_s);
            std::cerr << "node_iitree_ovlp " << idx << " " << start_in_s << " " << end_in_s << " " << pos_to_string(pos_start_in_q) << " " << pos_to_string(pos_end_in_q) << std::endl;
            // get the positions in s on either side of this range in Q
            // note the need to handle orientation
            bool is_rev_pos = is_rev(pos_start_in_q);
            std::vector<size_t> path_before_ovlp, path_after_ovlp;
            uint64_t start_in_q = (is_rev_pos ? offset(pos_end_in_q) : offset(pos_start_in_q));
            uint64_t end_in_q = (is_rev_pos ? offset(pos_start_in_q) : offset(pos_end_in_q));
            path_iitree.overlap(start_in_q-1, start_in_q, path_before_ovlp);
            path_iitree.overlap(end_in_q, end_in_q+1, path_after_ovlp);
            // and map these back into positions and orientations in S
            for (auto& idx : path_before_ovlp) {
                uint64_t start = path_iitree.start(idx);
                uint64_t end = path_iitree.end(idx);
                pos_t pos_start_in_s = path_iitree.data(idx);
                // if o_ur end is equal to start_in_q, it'd be a non-trivial edge
                std::cerr << "path_iitree before ovlp end == " << end << " ? " << start_in_q << " pos in s " << pos_to_string(pos_start_in_s) << std::endl;
            }
            for (auto& idx : path_after_ovlp) {
                uint64_t start = path_iitree.start(idx);
                uint64_t end = path_iitree.end(idx);
                pos_t pos_start_in_s = path_iitree.data(idx);
                // if o_ur end is equal to start_in_q, it'd be a non-trivial edge
                std::cerr << "path_iitree after ovlp end == " << end << " ? " << start_in_q << " pos in s " << pos_to_string(pos_start_in_s) << std::endl;
            }            
        }
        // and find their neighbors by looking in the path_iitree
        
        // decide and record which links these imply
    }
    link_mmset.index();
}

}
