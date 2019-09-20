#include "transclosure.hpp"
#include "spinlock.hpp"

namespace seqwish {

size_t compute_transitive_closures(
    seqindex_t& seqidx,
    mmmulti::iitree<uint64_t, pos_t>& aln_iitree, // input alignment matches between query seqs
    const std::string& seq_v_file,
    mmmulti::iitree<uint64_t, pos_t>& node_iitree, // maps graph seq ranges to input seq ranges
    mmmulti::iitree<uint64_t, pos_t>& path_iitree, // maps input seq ranges to graph seq ranges
    uint64_t repeat_max,
    uint64_t min_transclose_len) {
    // open seq_v_file
    std::ofstream seq_v_out(seq_v_file.c_str());
    // remember the elements of Q we've seen
    sdsl::bit_vector q_seen_bv(seqidx.seq_length());
    // for each base in seqidx
    //   collect mapped bases
    //   mark our seen bitvector
    //   emit the base
    //   record entries in node_mm and path_mm
    uint64_t input_seq_length = seqidx.seq_length();
    /*
    std::cerr << "all overlaps" << std::endl;
    std::vector<size_t> all_ovlp;
    aln_iitree.overlap(0, input_seq_length, all_ovlp);
    for (auto& s : all_ovlp) {
        uint64_t start = aln_iitree.start(s);
        uint64_t end = aln_iitree.end(s);
        pos_t pos = aln_iitree.data(s);
        std::cerr << "overlap " << start << "-" << end
                  << " position offset " << offset(pos) << (is_rev(pos)?"-":"+") << std::endl;
    }
    */
    // a buffer of ranges to write into our iitree, arranged by range ending position in Q
    // we flush those intervals that don't get extended into the next position in S
    //
    // this maps from a position in Q (our input seqs concatenated, offset and orientation)
    // to a range (start and length) in S (our graph sequence vector)
    std::map<pos_t, std::pair<uint64_t, uint64_t>> range_buffer;
    // here we try to find a growing range to extend
    auto extend_range = [&](const uint64_t& s_pos, const pos_t& q_pos) {
        // find a position in the map that we can add onto
        pos_t q_last_pos = q_pos;
        decr_pos(q_last_pos);
        auto f = range_buffer.find(q_last_pos);
        // if one doesn't exist, add the range
        if (f == range_buffer.end()) {
            range_buffer[q_pos] = std::make_pair(s_pos, 1);
        } else {
            // if one does, we expand the range and drop it back into the map at the new Q end pos
            std::pair<uint64_t, uint64_t> x = f->second;
            range_buffer.erase(f);
            ++x.second; // increment the match length
            range_buffer[q_pos] = x; // and stash it
        }
    };
    auto flush_ranges = [&](const uint64_t& s_pos) {
        // for each range, we're going to see if we've stepped more than one past the end
        // if we have, we'll write them out
        std::map<pos_t, std::pair<uint64_t, uint64_t>>::iterator it = range_buffer.begin();
        while (it != range_buffer.end()) {
            if (it->second.first + it->second.second != s_pos) {
                // flush
                std::cerr << "flushing " << offset(it->first) << ":" << (is_rev(it->first) ? "-" : "+")
                          << " -> " << it->second.first << ":" << it->second.second << std::endl;
                //node_mm.append(seq_v_length, j);
                //path_mm.append(offset(j), make_pos_t(seq_v_length,is_rev(j)));
                uint64_t match_length = it->second.second;
                uint64_t match_start_in_s = it->second.first;
                uint64_t match_end_in_s = match_start_in_s + match_length;
                pos_t match_pos_in_q = it->first;
                uint64_t match_start_in_q = offset(match_pos_in_q);
                uint64_t match_end_in_q = match_start_in_q + match_length;
                pos_t match_pos_in_s = make_pos_t(match_start_in_s, is_rev(match_pos_in_q));
                if (is_rev(match_pos_in_q)) std::swap(match_start_in_q, match_end_in_q);
                node_iitree.add(match_start_in_s, match_end_in_s, match_pos_in_q);
                path_iitree.add(match_start_in_q, match_end_in_q, match_pos_in_s);
                it = range_buffer.erase(it);
            } else {
                ++it;
            }
        }
    };
    //std::cerr << "input seq len " << input_seq_length << std::endl;
    // new way, accumulate pairs of ranges in S (output seq of graph) and Q (input seqs concatenated)
    // but how
    // maybe we set up a handler that decides when to flush things
    // we flush when we stop extending
    // we determine when we stop extending when we have stepped a bp and broke our range extension
    for (uint64_t i = 1; i <= input_seq_length; ++i) {
        //std::cerr << q_seen_bv << std::endl;
        if (q_seen_bv[i-1]) continue;
        // write base
        char base = seqidx.at(i-1);
        seq_v_out << seqidx.at(i-1);
        size_t seq_v_length = seq_v_out.tellp();
        uint64_t flushed = range_buffer.size();
        flush_ranges(seq_v_length);
        flushed -= range_buffer.size();
        std::cerr << "seq " << seq_v_length << " " << base << " buf size " << range_buffer.size() << " flushed " << flushed << std::endl;
        // mark current 
        q_seen_bv[i-1] = 1;
        // emit current
        std::set<std::pair<pos_t, uint64_t>> todo;
        std::unordered_map<uint64_t, uint64_t> seen_seqs;
        todo.insert(std::make_pair(make_pos_t(i, false), min_transclose_len));
        seen_seqs[seqidx.seq_id_at(offset(i))]++;
        while (!todo.empty()) {
            pos_t j = todo.begin()->first;
            uint64_t match_len = todo.begin()->second;
            todo.erase(todo.begin());
            //assert(q_seen_bv[offset(j)-1]==1);
            extend_range(seq_v_length, j);
            // todo rip this out
            if (min_transclose_len && match_len < min_transclose_len) {
                continue;
            }
            std::vector<size_t> ovlp;
            uint64_t n = offset(j);
            //std::cerr << "offset j " << n << std::endl;
            aln_iitree.overlap(n, n+1, ovlp);
            for (auto& s : ovlp) {
                uint64_t start = aln_iitree.start(s);
                uint64_t end = aln_iitree.end(s);
                pos_t pos = aln_iitree.data(s);
                //std::cerr << " with overlap " << start << "-" << end << std::endl;
                //std::cerr << "and position offset " << offset(pos) << (is_rev(pos)?"-":"+") << std::endl;
                //std::cerr << "n " << n << " start " << start << std::endl;
                incr_pos(pos, n - start);
                uint64_t k = offset(pos);
                if (k && !q_seen_bv[k-1]) {
                    uint64_t seq_id = seqidx.seq_id_at(offset(pos));
                    auto& c = seen_seqs[seq_id];
                    if (!repeat_max || c < repeat_max) {
                        ++c;
                        q_seen_bv[k-1] = 1;
                        todo.insert(std::make_pair(make_pos_t(offset(pos),is_rev(pos)^is_rev(j)), end - start));
                    }
                }
            }
        }
    }
    // close the graph sequence vector
    size_t seq_bytes = seq_v_out.tellp();
    seq_v_out.close();
    flush_ranges(seq_bytes+2);
    // build node_mm and path_mm indexes
    node_iitree.index();
    path_iitree.index();
    return seq_bytes;
}

}
