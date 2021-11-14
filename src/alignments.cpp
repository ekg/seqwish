#include "alignments.hpp"

namespace seqwish {

void paf_worker(
    igzstream& paf_in,
    std::atomic<bool>& paf_more,
    std::mutex& paf_in_mutex,
    mmmulti::iitree<uint64_t, pos_t>& aln_iitree,
    const seqindex_t& seqidx,
    const uint64_t& min_match_len) {
    while (paf_more.load()) {
        paf_in_mutex.lock();
        std::string line;
        std::getline(paf_in, line);
        paf_more.store(paf_in.good());
        paf_in_mutex.unlock();
        // Check if there is something to parse
        if (line.empty()) break;
        paf_row_t paf(line);
        // Check if the coordinates are reasonable
        if (paf.query_sequence_length == 0 || paf.target_sequence_length == 0 ||
            // Query/Target start (0-based; BED-like; closed)
            paf.query_start >= paf.query_sequence_length || paf.query_end > paf.query_sequence_length || paf.query_start >= paf.query_end ||
            // Query/Target end (0-based; BED-like; open)
            paf.target_start >= paf.target_sequence_length || paf.target_end > paf.target_sequence_length || paf.target_start >= paf.target_end) break;
        size_t query_idx = seqidx.rank_of_seq_named(paf.query_sequence_name);
        //size_t query_len = seqidx.nth_seq_length(query_idx);
        size_t target_idx = seqidx.rank_of_seq_named(paf.target_sequence_name);
        //size_t target_len = seqidx.nth_seq_length(target_idx);
        bool q_rev = !paf.query_target_same_strand;
        size_t q_all_pos = (q_rev ? seqidx.pos_in_all_seqs(query_idx, paf.query_end, false) - 1
                            : seqidx.pos_in_all_seqs(query_idx, paf.query_start, false));
        size_t t_all_pos = seqidx.pos_in_all_seqs(target_idx, paf.target_start, false);
        pos_t q_pos = make_pos_t(q_all_pos, q_rev);
        pos_t t_pos = make_pos_t(t_all_pos, false);
        for (auto& c : paf.cigar) {
            switch (c.op) {
            case 'M':
            case '=':
            case 'X':
            {
                pos_t q_pos_match_start = q_pos;
                pos_t t_pos_match_start = t_pos;
                uint64_t match_len = 0;
                auto add_match =
                    [&](void) {
                        if (match_len && match_len >= min_match_len) {
                            if (is_rev(q_pos)) {
                                pos_t x_pos = q_pos;
                                decr_pos(x_pos); // to guard against underflow when our start is 0-, we need to decr in pos_t space
                                aln_iitree.add(offset(x_pos), offset(q_pos_match_start)+1, make_pos_t(offset(t_pos)-1, true));
                                aln_iitree.add(offset(t_pos_match_start), offset(t_pos), make_pos_t(offset(q_pos_match_start), true));
                            } else {
                                aln_iitree.add(offset(q_pos_match_start), offset(q_pos), t_pos_match_start);
                                aln_iitree.add(offset(t_pos_match_start), offset(t_pos), q_pos_match_start);
                            }
                        }
                    };
                for (size_t i = 0; i < c.len; ++i) {
                    char query_base = seqidx.at_pos(q_pos);
                    char target_base = seqidx.at_pos(t_pos);
                    if (query_base == target_base
                        && query_base != 'N'
                        && offset(q_pos) != offset(t_pos)) { // guard against self mappings
                        if (match_len == 0) {
                            q_pos_match_start = q_pos;
                            t_pos_match_start = t_pos;
                        }
                        ++match_len;
                        incr_pos(q_pos);
                        incr_pos(t_pos);
                    } else {
                        add_match();
                        incr_pos(q_pos);
                        incr_pos(t_pos);
                        match_len = 0;
                        // break out the last match
                    }
                }
                // handle any last match
                add_match();
            }
                break;
            case 'I':
                //std::cerr << "ins " << c.len << std::endl;
                incr_pos(q_pos, c.len);
                break;
            case 'D':
                //std::cerr << "del " << c.len << std::endl;
                incr_pos(t_pos, c.len);
                break;
            default: break;
            }
        }
    }
}


void unpack_paf_alignments(const std::string& paf_file,
                           mmmulti::iitree<uint64_t, pos_t>& aln_iitree,
                           const seqindex_t& seqidx,
                           const uint64_t& min_match_len,
                           const uint64_t& num_threads) {
    // go through the PAF file
    igzstream paf_in(paf_file.c_str());
    if (!paf_in.good()) {
        std::cerr << "[seqwish::alignments] error: PAF file " << paf_file << " is not good!" << std::endl;
        exit(1);
    }
    std::mutex paf_in_mutex;
    std::atomic<bool> paf_more; paf_more.store(true);
    std::vector<std::thread> workers; workers.reserve(num_threads);
    for (uint64_t t = 0; t < num_threads; ++t) {
        workers.emplace_back(paf_worker, std::ref(paf_in), std::ref(paf_more), std::ref(paf_in_mutex), std::ref(aln_iitree), std::ref(seqidx), std::ref(min_match_len));
    }
    for (uint64_t t = 0; t < num_threads; ++t) {
        workers[t].join();
    }
}

}
