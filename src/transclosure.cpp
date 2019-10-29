#include "transclosure.hpp"
#include "spinlock.hpp"

namespace seqwish {

size_t compute_transitive_closures(
    seqindex_t& seqidx,
    range_pos_iitii& aln_iitree, // input alignment matches between query seqs
    const std::string& seq_v_file,
    range_pos_iitii::builder& node_iitree_builder,
    range_pos_iitii::builder& path_iitree_builder,
    uint64_t repeat_max,
    uint64_t min_transclose_len) {
    // open seq_v_file
    std::ofstream seq_v_out(seq_v_file.c_str());
    // remember the elements of Q we've seen
    sdsl::bit_vector q_seen_bv(seqidx.seq_length());
    uint64_t input_seq_length = seqidx.seq_length();
    // a buffer of ranges to write into our iitree, arranged by range ending position in Q
    // we flush those intervals that don't get extended into the next position in S
    // this maps from a position in Q (our input seqs concatenated, offset and orientation)
    // to a range (start and length) in S (our graph sequence vector)
    // we are mapping from the /last/ position in the matched range, not the first
    std::map<pos_t, std::pair<uint64_t, uint64_t>> range_buffer;
    // here we try to find a growing range to extend
    auto extend_range = [&](const uint64_t& s_pos, const pos_t& q_pos) {
        // find a position in the map that we can add onto
        // it must match position and orientation
        pos_t q_last_pos = q_pos;
        decr_pos(q_last_pos);
        auto f = range_buffer.find(q_last_pos);
        // if one doesn't exist, add the range
        if (f == range_buffer.end()) {
            range_buffer[q_pos] = std::make_pair(s_pos, 1);
        } else {
            // if one does, check that it matches our extension,
            std::pair<uint64_t, uint64_t> x = f->second;
            if (x.first + x.second == s_pos) {
                // if so we expand its range and drop it back into the map at the new Q end pos
                range_buffer.erase(f);
                ++x.second; // increment the match length
                range_buffer[q_pos] = x; // and stash it
            } else {
                // if it doesn't, we store a new range
                range_buffer[q_pos] = std::make_pair(s_pos, 1);
            }
        }
    };
    auto flush_ranges = [&](const uint64_t& s_pos) {
        // for each range, we're going to see if we've stepped more than one past the end
        // if we have, we'll write them out
        std::map<pos_t, std::pair<uint64_t, uint64_t>>::iterator it = range_buffer.begin();
        while (it != range_buffer.end()) {
            if (it->second.first + it->second.second != s_pos) {
                // flush
                uint64_t match_length = it->second.second;
                uint64_t match_start_in_s = it->second.first;
                uint64_t match_end_in_s = match_start_in_s + match_length;
                pos_t match_end_pos_in_q = it->first;
                bool is_rev_match = is_rev(match_end_pos_in_q);
                // TODO appreciate why you're doing this here, presumably it's because of how we're running the transclosure
                if (!is_rev_match) incr_pos(match_end_pos_in_q, 1);
                //incr_pos(match_end_pos_in_q, 1);
                pos_t match_start_pos_in_q = match_end_pos_in_q;
                decr_pos(match_start_pos_in_q, match_length);
                uint64_t match_end_in_q = offset(match_end_pos_in_q);
                uint64_t match_start_in_q = offset(match_start_pos_in_q);
                pos_t match_pos_in_s = make_pos_t(match_start_in_s, is_rev_match);
                pos_t match_pos_in_q = make_pos_t(match_start_in_q, is_rev_match);
                if (is_rev_match) {
                    // go form transclosure model to the same pattern we have in the alignment iitree
                    // 1-based half open intervals, positions map to start and orientation in S and Q
                    std::swap(match_start_in_q, match_end_in_q);
                    incr_pos(match_pos_in_q, 1);
                    decr_pos(match_pos_in_s, match_length - 1);
                }
                node_iitree_builder.add({match_start_in_s, match_end_in_s, match_pos_in_q});
                path_iitree_builder.add({match_start_in_q, match_end_in_q, match_pos_in_s});
                it = range_buffer.erase(it);
            } else {
                ++it;
            }
        }
    };
    // range compressed model
    // accumulate pairs of ranges in S (output seq of graph) and Q (input seqs concatenated)
    // we flush when we stop extending
    // we determine when we stop extending when we have stepped a bp and broke our range extension
    uint64_t last_seq_id = seqidx.seq_id_at(1);
    for (uint64_t i = 1; i <= input_seq_length; ++i) {
        if (q_seen_bv[i-1]) continue;
        // write base
        char base = seqidx.at(i-1);
        seq_v_out << seqidx.at(i-1);
        size_t seq_v_length = seq_v_out.tellp();
        //uint64_t flushed = range_buffer.size();
        uint64_t curr_seq_id = seqidx.seq_id_at(i);
        if (curr_seq_id != last_seq_id) {
            flush_ranges(seq_v_length+1); // hack to force flush at sequence change
            last_seq_id = curr_seq_id;
        } else {
            flush_ranges(seq_v_length);
        }
        //flushed -= range_buffer.size();
        //std::cerr << "seq " << seq_v_length << " " << base << " buf size " << range_buffer.size() << " flushed " << flushed << std::endl;
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
            // optionally require a minimum length of match to transclose through
            if (min_transclose_len && match_len < min_transclose_len) {
                continue;
            }
            uint64_t n = offset(j);
            //std::cerr << "offset j " << n << std::endl;
            std::vector<range_pos_t> ovlp = aln_iitree.overlap(n, n+1);
            for (auto& s : ovlp) {
                auto& start = s.start;
                auto& end = s.end;
                auto& pos = s.pos;
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
    assert(range_buffer.empty());
    return seq_bytes;
}

}
