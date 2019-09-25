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
    for (size_t i = 0; i < graph_size; ++i) {
        size_t j = i+1;
        uint64_t from = i;
        uint64_t to = j;
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
    size_t num_seqs = seqidx.n_seqs();
    for (size_t i = 1; i <= num_seqs; ++i) {
        std::cerr << "at seq " << i << std::endl;
        size_t j = seqidx.nth_seq_offset(i)+1;
        size_t k = j+seqidx.nth_seq_length(i)-1;
        std::cerr << "j " << j << " k " << k << std::endl;
        std::vector<size_t> ovlp_start;
        path_iitree.overlap(j, j+1, ovlp_start);
        assert(ovlp_start.size() == 1);
        pos_t p1 = path_iitree.data(ovlp_start.front());
        if (is_rev(p1)) {
            seq_id_bv[offset(p1)] = 1;
        } else {
            std::cerr << "offset(p1) = " << offset(p1) << std::endl;
            seq_id_bv[offset(p1)-1] = 1;
        }
        std::vector<size_t> ovlp_end;
        path_iitree.overlap(k, k+1, ovlp_end);
        assert(ovlp_end.size() == 1);
        pos_t p2 = path_iitree.data(ovlp_end.front());
        if (is_rev(p2)) {
            seq_id_bv[offset(p2)-1] = 1;
        } else {
            std::cerr << "offset(p2) = " << offset(p2) << std::endl;
            seq_id_bv[offset(p2)] = 1;
        }
    }
}

}
