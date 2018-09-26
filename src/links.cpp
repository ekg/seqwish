#include "links.hpp"
#include "dmultimap.hpp"

namespace seqwish {

void derive_links(seqindex_t& seqidx,
                  size_t graph_length,
                  dmultimap<uint64_t, pos_t>& path_mm,
                  dmultimap<pos_t, pos_t>& link_fwd_mm,
                  dmultimap<pos_t, pos_t>& link_rev_mm) {
    // rewrite the sequences in seqidx as pairs of positions in S
    // for each sequence in seqidx
    size_t num_seqs = seqidx.n_seqs();
#pragma omp parallel for
    for (size_t i = 1; i <= num_seqs; ++i) {
        size_t j = seqidx.nth_seq_offset(i);
        size_t k = j+seqidx.nth_seq_length(i);
        //std::cerr << seqidx.nth_name(i) << " " << seqidx.nth_seq_length(i) << " " << j << " " << k << std::endl;
        for ( ; j < k-1; ++j) {
            std::vector<pos_t> v1 = path_mm.unique_values(j+1);
            std::vector<pos_t> v2 = path_mm.unique_values(j+2);
            // each input base should only map one place in the graph
            assert(v1.size() == v2.size() == 1);
            auto& p1 = v1.front();
            auto& p2 = v2.front();
            //std::cerr << pos_to_string(p1) << " " << pos_to_string(p2) << std::endl;
            link_fwd_mm.append(p1, p2);
            link_fwd_mm.append(rev_pos_t(p2), rev_pos_t(p1));
            link_rev_mm.append(p2, p1);
            link_rev_mm.append(rev_pos_t(p1), rev_pos_t(p2));
        }
    }
    link_fwd_mm.index(make_pos_t(graph_length, true));
    link_rev_mm.index(make_pos_t(graph_length, true));
}

}
