#include "compact.hpp"
#include "mmmultimap.hpp"

namespace seqwish {

void compact_nodes(
    seqindex_t& seqidx,
    size_t graph_size,
    mmmulti::iitree<uint64_t, pos_t>& node_iitree,
    mmmulti::iitree<uint64_t, pos_t>& path_iitree,
    sdsl::bit_vector& seq_id_bv) {
    //seq_id_bv;
    // for each pair of positions in the graph base seq
    // do we have any links from the first that don't go to the second?
    // do we have any links to the second that don't come from the first?
    seq_id_bv[0] = 1; // set first node start
#pragma omp parallel for schedule(dynamic)
    for (size_t i = 1; i <= graph_size; ++i) {
        uint64_t from = i;
        uint64_t to = i + 1;
        std::vector<size_t> from_ovlp, to_ovlp;
        node_iitree.overlap(from, from+1, from_ovlp);
        node_iitree.overlap(to, to+1, to_ovlp);
        // query intervals overlapping this position in the node_iitree
        std::cerr << "in compact " << from << " -> " << to << std::endl;
        std::cerr << from_ovlp.size() << " " << to_ovlp.size() << std::endl;
        if (from_ovlp ==  to_ovlp) {
            std::cerr << "continuing" << std::endl;
        } else {
            // mark a node start
#pragma omp critical (seq_id_bv)
            {
                std::cerr << "marking node start " << i << " of " << seq_id_bv.size() << std::endl;
                seq_id_bv[i] = 1;
            }
        }
    }
    std::cerr << graph_size << " " << seq_id_bv.size() << std::endl;
    seq_id_bv[graph_size] = 1;
}

}
