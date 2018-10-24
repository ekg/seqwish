#include "alignments.hpp"

namespace seqwish {

void unpack_alignments(const std::string& paf_file,
                       dmultimap<uint64_t, pos_t>& aln_mm,
                       seqindex_t& seqidx) {
    // go through the PAF file
    
    igzstream paf_in(paf_file.c_str());
    if (!paf_in.good()) assert("PAF is not good!");
    std::string line;
    std::vector<std::string> lines;
    while (std::getline(paf_in, line)) {
        lines.push_back(line);
    }
#pragma omp parallel for schedule(dynamic)
    for (size_t i = 0; i < lines.size(); ++i) {
        auto& line = lines[i];
        paf_row_t paf(line);
        //std::cerr << paf << std::endl;
        size_t query_idx = seqidx.rank_of_seq_named(paf.query_sequence_name);
        size_t query_len = seqidx.nth_seq_length(query_idx);
        size_t target_idx = seqidx.rank_of_seq_named(paf.target_sequence_name);
        size_t target_len = seqidx.nth_seq_length(target_idx);
        bool q_rev = !paf.query_target_same_strand;
        //std::cerr << "query_idx " << query_idx << " " << seqidx.nth_seq_length(query_idx) << std::endl;
        //std::cerr << "target_idx " << target_idx << " " << seqidx.nth_seq_length(target_idx) << std::endl;
        size_t target_start = paf.target_start;
        size_t query_start = q_rev ? seqidx.nth_seq_length(query_idx)-paf.query_end : paf.query_start;
        //std::cerr << query_start << " " << target_start << std::endl;

        // these calls convert to 1-based positions as 0 has a special meaning in the dmultimap 
        size_t q_all_pos = 1 + seqidx.pos_in_all_seqs(query_idx, query_start, q_rev);
        size_t t_all_pos = 1 + seqidx.pos_in_all_seqs(target_idx, target_start, false);
        //std::cerr << "q_all_pos " << q_all_pos << std::endl;
        //std::cerr << "t_all_pos " << t_all_pos << std::endl;
        pos_t q_pos = make_pos_t(q_all_pos, q_rev);
        pos_t t_pos = make_pos_t(t_all_pos, false);
        for (auto& c : paf.cigar) {
            switch (c.op) {
            case 'M':
                //std::cerr << "match " << c.len << std::endl;
                for (size_t i = 0; i < c.len; ++i) {
                    /*
                    std::cerr << seqidx.at_pos(t_pos) << " vs " << seqidx.at_pos(q_pos) << " ... "
                              << offset(t_pos) << " " << offset(q_pos)
                              << std::endl;
                    */
                    if (seqidx.at_pos(q_pos) == seqidx.at_pos(t_pos)
                        && offset(q_pos) != offset(t_pos)) {
                        assert(offset(q_pos) && offset(t_pos));
                        //aln_mm.append(offset(q_pos), { make_pos_t(offset(t_pos), q_rev), paf.bases_in_mapping });
                        //aln_mm.append(offset(t_pos), { make_pos_t(offset(q_pos), q_rev), paf.bases_in_mapping });
                        aln_mm.append(offset(q_pos), make_pos_t(offset(t_pos), q_rev));
                        aln_mm.append(offset(t_pos), make_pos_t(offset(q_pos), q_rev));
                        //} else {
                        //std::cerr << "mismatch" << std::endl;
                    }
                    incr_pos(q_pos);
                    incr_pos(t_pos);
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
    aln_mm.index(seqidx.seq_length());
}

/*
void filter_alignments(dmultimap<pos_t, aln_pos_t>& aln_mm,
                       dmultimap<pos_t, pos_t>& aln_filt_mm,
                       uint64_t aln_min_length,
                       uint64_t aln_keep_n_longest,
                       seqindex_t& seqidx) {
    for (size_t i = 1; i <= seqidx.seq_length(); ++i) {
        std::vector<aln_pos_t> at_pos = aln_mm.unique_values(i);
        std::sort(at_pos.begin(), at_pos.end(), [&](const aln_pos_t& a, const aln_pos_t& b){ return a.aln_length > b.aln_length; });
        uint64_t to_keep = (aln_keep_n_longest ? std::min(aln_keep_n_longest, (uint64_t)at_pos.size()) : at_pos.size());
        //std::cerr << "keep at " << i << ": ";
        //for (auto& k : at_pos) std::cerr << " " << k.aln_length;  std::cerr << std::endl;
        for (size_t j = 0; j < to_keep; ++j) {
            auto& a = at_pos[j];
            if (a.aln_length < aln_min_length) break;
            aln_filt_mm.append(i, at_pos[j].pos);
        }
    }
    aln_filt_mm.index(seqidx.seq_length());
}
*/

}
