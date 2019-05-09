#include "compact.hpp"
#include "mmmultimap.hpp"

namespace seqwish {

void compact_nodes(
    size_t graph_size,
    multimap<pos_t, pos_t>& link_fwd_mm,
    multimap<pos_t, pos_t>& link_rev_mm,
    sdsl::bit_vector& seq_id_bv) {
    //seq_id_bv;
    // for each pair of positions in the graph base seq
    // do we have any links from the first that don't go to the second?
    // do we have any links to the second that don't come from the first?
    seq_id_bv[0] = 1; // set first node start
#pragma omp parallel for schedule(dynamic)
    for (size_t i = 1; i < graph_size; ++i) {
        size_t j = i+1;
        pos_t from = make_pos_t(i, false);
        pos_t to = make_pos_t(j, false);
        std::vector<pos_t> from_first = link_fwd_mm.unique_values(from);
        std::vector<pos_t> to_second = link_rev_mm.unique_values(to);
        //std::cerr << "in compact " << pos_to_string(from) << " -> " << pos_to_string(to) << std::endl;
        //std::cerr << from_first.size() << " " << to_second.size() << std::endl;
        if (from_first.size() == 1 && from_first.front() == to
            && to_second.size() == 1 && to_second.front() == from) {
        } else {
            // mark a node start
#pragma omp critical (seq_id_bv)
            seq_id_bv[i] = 1;
        }
    }
}

}
