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
    //std::cerr << seq_id_cbv << std::endl;
    /*
    for (size_t i = 0; i < seq_id_cbv.size(); ++i) {
        std::cerr << i << " rank " << seq_id_cbv_rank(i) << std::endl;
    }
    */
#pragma omp parallel for schedule(dynamic)
    for (size_t id = 1; id <= n_nodes; ++id) {
        uint64_t node_start_in_s = seq_id_cbv_select(id); // node_iitree and path_iitree are 1-based
        uint64_t node_end_in_s = seq_id_cbv_select(id+1);
        //std::cerr << "links for node " << id << " start " << node_start_in_s << " end " << node_end_in_s << std::endl;
        // find the things on both sides of our node by looking in the node_iitree, finding what bits of the paths (in Q)
        // are there, and seeing what's on either side of them to decide what links we need
        std::vector<size_t> node_ovlp;
        node_iitree.overlap(node_start_in_s, node_end_in_s, node_ovlp);
        for (auto& idx : node_ovlp) {
            uint64_t ovlp_start_in_s = node_iitree.start(idx);
            uint64_t ovlp_end_in_s = node_iitree.end(idx);
            uint64_t ovlp_length = ovlp_end_in_s - ovlp_start_in_s;
            pos_t pos_start_in_q = node_iitree.data(idx);
            pos_t pos_end_in_q = pos_start_in_q;
            incr_pos(pos_end_in_q, ovlp_length - (ovlp_end_in_s - node_end_in_s));
            incr_pos(pos_start_in_q, node_start_in_s - ovlp_start_in_s);
            //std::cerr << "node_iitree_ovlp " << idx << " " << ovlp_start_in_s << " " << ovlp_end_in_s << " " << pos_to_string(pos_start_in_q) << " " << pos_to_string(pos_end_in_q) << std::endl;
            // get the positions in s on either side of this range in Q
            // note the need to handle orientation
            // determine which sequence we're in
            // find its boundaries
            bool curr_step_is_rev = is_rev(pos_start_in_q);
            uint64_t start_in_q = (curr_step_is_rev ? offset(pos_end_in_q)+1 : offset(pos_start_in_q));
            uint64_t end_in_q = (curr_step_is_rev ? offset(pos_start_in_q)+1 : offset(pos_end_in_q));
            //std::cerr << "start_in_q " << start_in_q << " end_in_q " << end_in_q << std::endl;
            uint64_t seq_id = seqidx.seq_id_at(start_in_q);
            uint64_t id_end = seqidx.seq_id_at(end_in_q-1);
            //std::cerr << "sequence id at start " << seq_id << " and end " << id_end << std::endl;
            assert(seq_id == id_end);
            
            uint64_t seq_start = seqidx.nth_seq_offset(seq_id);
            uint64_t seq_end = seq_start + seqidx.nth_seq_length(seq_id);
            //std::cerr << "seq boundaries " << seq_start << "-" << seq_end << std::endl;
            std::vector<size_t> path_before_ovlp, path_after_ovlp;
            // and only consider cases where we'd be within the boundaries
            if (start_in_q-1 >= seq_start) {
                path_iitree.overlap(start_in_q-1, start_in_q, path_before_ovlp);
                // map these back into positions and orientations in S
                for (auto& idx : path_before_ovlp) {
                    // this is the larger range we're in
                    uint64_t ovlp_start_in_q = path_iitree.start(idx);
                    uint64_t ovlp_end_in_q = path_iitree.end(idx);
                    pos_t pos_start_in_s = path_iitree.data(idx);
                    pos_t pos_end_in_s = pos_start_in_s;
                    if (ovlp_end_in_q > start_in_q) {
                        //std::cerr << "contained in previous range on q" << std::endl;
                        incr_pos(pos_end_in_s, (ovlp_end_in_q - ovlp_start_in_q) - (ovlp_end_in_q - start_in_q) - 1);
                    } else {
                        //std::cerr << "not contained in previous range on q" << std::endl;
                        incr_pos(pos_end_in_s, ovlp_end_in_q - ovlp_start_in_q - 1);
                    }
                    // now calculate where we should land with respect to our current node
                    // start_in_s and end_in_s be respecting
                    /*
                    std::cerr << "path_iitree before ovlp " << ovlp_start_in_q << ".." << ovlp_end_in_q << " vs " << start_in_q << ".." << end_in_q
                              << " prev pos_end_in_s " << pos_to_string(pos_end_in_s)
                              << " node boundaries " << node_start_in_s << ".." << node_end_in_s << std::endl;
                    */
                    // find which node we're in here, and record a link
                    uint64_t prev_id = seq_id_cbv_rank(offset(pos_end_in_s)+1);
                    //std::cerr << "id? " << prev_id << std::endl;
                    // relative orientations please
                    bool prev_step_is_rev = is_rev(pos_end_in_s);
                    //std::cerr << "link " << (prev_step_is_rev ? "-" : "+") << " " << (curr_step_is_rev ? "-" : "+") << std::endl;
                    /*
                    std::cerr << "link " << pos_to_string(make_pos_t(prev_id, prev_step_is_rev))
                              << " " << pos_to_string(make_pos_t(id, curr_step_is_rev)) << std::endl;
                    */
                    link_mmset.append(std::make_pair(make_pos_t(prev_id, prev_step_is_rev), make_pos_t(id, curr_step_is_rev)));
                }
            }
            if (end_in_q+1 <= seq_end) {
                path_iitree.overlap(end_in_q, end_in_q+1, path_after_ovlp);
                for (auto& idx : path_after_ovlp) {
                    uint64_t ovlp_start_in_q = path_iitree.start(idx);
                    uint64_t ovlp_end_in_q = path_iitree.end(idx);
                    pos_t pos_start_in_s = path_iitree.data(idx);
                    if (ovlp_start_in_q < end_in_q) {
                        //std::cerr << "contained in previous range on q" << std::endl;
                        incr_pos(pos_start_in_s, end_in_q - ovlp_start_in_q);
                        //} else {
                        //std::cerr << "not contained in previous range on q" << std::endl;
                    }
                    //incr_pos(pos_start_in_s, ovlp_end_in_q - ovlp_start_in_q);
                    /*
                    std::cerr << "path_iitree after ovlp " << ovlp_start_in_q << ".." << ovlp_end_in_q << " vs " << start_in_q << ".." << end_in_q
                              << " next pos_start_in_s " << pos_to_string(pos_start_in_s)
                              << " node boundaries " << node_start_in_s << ".." << node_end_in_s << std::endl;
                    */
                    // find which node we're in here, and record a link
                    uint64_t next_id = seq_id_cbv_rank(offset(pos_start_in_s)+1);
                    //std::cerr << "id? " << next_id << std::endl;
                    bool next_step_is_rev = is_rev(pos_start_in_s);
                    //std::cerr << "link " << (curr_step_is_rev ? "-" : "+") << " " << (next_step_is_rev ? "-" : "+") << std::endl;
                    /*
                    std::cerr << "link " << pos_to_string(make_pos_t(id, curr_step_is_rev))
                              << " " << pos_to_string(make_pos_t(next_id, next_step_is_rev)) << std::endl;
                    */
                    link_mmset.append(std::make_pair(make_pos_t(id, curr_step_is_rev), make_pos_t(next_id, next_step_is_rev)));
                }
            }
        }
        // and find their neighbors by looking in the path_iitree
        
        // decide and record which links these imply
    }
    link_mmset.index();
}

}
