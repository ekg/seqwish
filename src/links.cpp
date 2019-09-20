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
        size_t node_start = seq_id_cbv_select(id);
        size_t node_end = seq_id_cbv_select(id+1);
        // find the things on both sides of our node by looking in 
        // link them together using their index in Q
        // if they are one bp apart in Q then by definition they should induce a link in the graph
    }
    link_mmset.index();
}

}
