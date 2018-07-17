#include "alignments.hpp"

namespace seqwish {

void unpack_alignments(const std::string& paf_file,
                       dmultimap<uint64_t, pos_t>& aln_mm,
                       seqindex_t& seqidx) {
    // go through the PAF file
    std::ifstream paf_in(paf_file.c_str());
    std::string line;
    while (std::getline(paf_in, line)) {
        paf_row_t paf(line);
        //std::cerr << paf << std::endl;
        size_t query_idx = seqidx.rank_of_seq_named(paf.query_sequence_name);
        size_t query_len = seqidx.nth_seq_length(query_idx);
        size_t target_idx = seqidx.rank_of_seq_named(paf.target_sequence_name);
        size_t target_len = seqidx.nth_seq_length(target_idx);
        bool t_rev = !paf.query_target_same_strand;
        //std::cerr << "query_idx " << query_idx << " " << seqidx.nth_seq_length(query_idx) << std::endl;
        //std::cerr << "target_idx " << target_idx << " " << seqidx.nth_seq_length(target_idx) << std::endl;
        // these calls convert to 1-based positions as 0 has a special meaning in the dmultimap 
        size_t q_all_pos = 1 + seqidx.pos_in_all_seqs(query_idx, paf.query_start, false);
        size_t t_all_pos = 1 + seqidx.pos_in_all_seqs(target_idx, paf.target_start, t_rev) - (t_rev?1:0);
        //std::cerr << "q_all_pos " << q_all_pos << std::endl;
        //std::cerr << "t_all_pos " << t_all_pos << std::endl;
        uint64_t q_pos = q_all_pos;
        pos_t t_pos = make_pos_t(t_all_pos, t_rev);
        for (auto& c : paf.cigar) {
            switch (c.op) {
            case 'M':
                for (size_t i = 0; i < c.len; ++i) {
                    //std::cerr << q_pos << ":" << pos_to_string(t_pos) << "/" << t_pos << std::endl;
                    aln_mm.append(q_pos, t_pos);
                    ++q_pos;
                    incr_pos(t_pos);
                }
                break;
            case 'I':
                q_pos += c.len;
                break;
            case 'D':
                incr_pos(t_pos, c.len);
                break;
            default: break;
            }
        }
    }
    aln_mm.index(seqidx.seq_length());
}

}
