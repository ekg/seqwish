#include "transclosure.hpp"
#include "mmmultimap.hpp"
#include "mmiitree.hpp"
#include "spinlock.hpp"

namespace seqwish {

size_t compute_transitive_closures(
    seqindex_t& seqidx,
    mmmulti::iitree<uint64_t, std::pair<pos_t, uint64_t>>& aln_iitree,
    const std::string& seq_v_file,
    mmmulti::map<uint64_t, pos_t>& node_mm,
    mmmulti::map<uint64_t, pos_t>& path_mm,
    uint64_t repeat_max,
    uint64_t min_transclose_len) {
    // open seq_v_file
    std::ofstream seq_v_out(seq_v_file.c_str());
    // remember the elements of Q we've seen
    sdsl::bit_vector q_seen_bv(seqidx.seq_length());
    //size_t seq_bytes = 0;
    // for each base in seqidx
    //   collect mapped bases
    //   mark our seen bitvector
    //   emit the base
    //   record entries in node_mm and path_mm
    uint64_t input_seq_length = seqidx.seq_length();
    for (uint64_t i = 1; i <= input_seq_length; ++i) {
        //std::cerr << "top " << i << std::endl;
        if (q_seen_bv[i-1]) continue;
        // write base
        //std::cerr << "writing seq out " << seqidx.at(i-1) << std::endl;
        seq_v_out << seqidx.at(i-1);
        size_t seq_v_length = seq_v_out.tellp();
        // mark current
        q_seen_bv[i-1] = 1;
        // emit current
        //node_mm.append(i, make_pos_t(i,false));
        //path_mm.append(i, make_pos_t(i,false));
        //std::cerr << "closure " << i << " to " << seq_v_length << std::endl;
        std::set<std::pair<pos_t, uint64_t>> todo;
        std::unordered_map<uint64_t, uint64_t> seen_seqs;
        todo.insert(std::make_pair(make_pos_t(i, false), min_transclose_len));
        seen_seqs[seqidx.seq_id_at(offset(i))]++;
        while (!todo.empty()) {
            pos_t j = todo.begin()->first;
            uint64_t match_len = todo.begin()->second;
            /*
            std::cerr << seq_v_length << "\t"
                      << offset(j) << "\t"
                      << (is_rev(j)?"-":"+") << std::endl;
            */
            todo.erase(todo.begin());
            //std::cerr << todo.size() << std::endl;
            //assert(q_seen_bv[offset(j)-1]==1);
            node_mm.append(seq_v_length, j);
            path_mm.append(offset(j), make_pos_t(seq_v_length,is_rev(j)));
            // todo rip this out
            if (min_transclose_len && match_len < min_transclose_len) {
                continue;
            }
            std::vector<size_t> ovlp;
            uint64_t n = offset(j);
            aln_iitree.overlap(n, n+1, ovlp);
            for (auto& s : ovlp) {
                uint64_t start = aln_iitree.start(s);
                uint64_t end = aln_iitree.end(s);
                const std::pair<pos_t, uint64_t>& interval = aln_iitree.data(s);
                // if it's long enough, include it
                // todo, use paramater for local smoothing
                // find the position in the match
                pos_t pos = interval.first;
                incr_pos(pos, n - start);
                uint64_t k = offset(pos);
                if (k && !q_seen_bv[k-1]) {
                    //std::cerr << "closing " << k << std::endl;
                    uint64_t seq_id = seqidx.seq_id_at(offset(pos));
                    //std::cerr << "seq id " << seq_id << std::endl;
                    auto& c = seen_seqs[seq_id];
                    if ((!repeat_max || c < repeat_max)) {
                        ++c;
                        q_seen_bv[k-1] = 1;
                        todo.insert(std::make_pair(make_pos_t(offset(pos),is_rev(pos)^is_rev(j)), end - start));
                    }
                }
            }
        }
        /*
        for (auto& c : seen_seqs) {
            std::cerr << "seen " << c.first << " " << c.second << std::endl;
        }
        */
    }
    // close the graph sequence vector
    size_t seq_bytes = seq_v_out.tellp();
    seq_v_out.close();
    // build node_mm and path_mm indexes
    node_mm.index(seq_bytes);
    path_mm.index(input_seq_length);
    return seq_bytes;
}

}
