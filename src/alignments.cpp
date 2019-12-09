#include "alignments.hpp"

namespace seqwish {

void unpack_paf_alignments(const std::string& paf_file,
                           range_pos_iitii::builder& aln_iitree_builder,
                           seqindex_t& seqidx,
                           uint64_t min_match_len) {
    // go through the PAF file
    igzstream paf_in(paf_file.c_str());
    if (!paf_in.good()) assert("PAF is not good!");
    uint64_t lines = std::count(std::istreambuf_iterator<char>(paf_in), 
                                std::istreambuf_iterator<char>(), '\n');
    paf_in.close();
    paf_in.open(paf_file.c_str());
//#pragma omp parallel for schedule(dynamic) // why is this broken now?
    for (size_t i = 0; i < lines; ++i) {
        std::string line;
#pragma omp critical (paf_in)
        std::getline(paf_in, line);
        paf_row_t paf(line);
        size_t query_idx = seqidx.rank_of_seq_named(paf.query_sequence_name);
        size_t query_len = seqidx.nth_seq_length(query_idx);
        size_t target_idx = seqidx.rank_of_seq_named(paf.target_sequence_name);
        size_t target_len = seqidx.nth_seq_length(target_idx);
        bool q_rev = !paf.query_target_same_strand;
        size_t q_all_pos = (q_rev ? seqidx.pos_in_all_seqs(query_idx, paf.query_end, false) - 1
                            : seqidx.pos_in_all_seqs(query_idx, paf.query_start, false));
        size_t t_all_pos = seqidx.pos_in_all_seqs(target_idx, paf.target_start, false);
        pos_t q_pos = make_pos_t(q_all_pos, q_rev);
        pos_t t_pos = make_pos_t(t_all_pos, false);
        for (auto& c : paf.cigar) {
            switch (c.op) {
            case 'M':
            {
                pos_t q_pos_match_start = q_pos;
                pos_t t_pos_match_start = t_pos;
                uint64_t match_len = 0;
                auto add_match = [&](void) {
                    if (match_len && match_len >= min_match_len) {
                        if (is_rev(q_pos)) {
                            aln_iitree_builder.add({offset(q_pos)+1, offset(q_pos_match_start)+1, make_pos_t(offset(t_pos)-1, true)});
                            aln_iitree_builder.add({offset(t_pos_match_start), offset(t_pos), make_pos_t(offset(q_pos_match_start), true)});
                        } else {
                            aln_iitree_builder.add({offset(q_pos_match_start), offset(q_pos), t_pos_match_start});
                            aln_iitree_builder.add({offset(t_pos_match_start), offset(t_pos), q_pos_match_start});
                        }
                    }
                };
                for (size_t i = 0; i < c.len; ++i) {
                    if (seqidx.at_pos(q_pos) == seqidx.at_pos(t_pos)
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
    //aln_iitree.index();
}

void unpack_sxs_alignments(const std::string& sxs_file,
                           range_pos_iitii::builder& aln_iitree_builder,
                           seqindex_t& seqidx,
                           uint64_t min_match_len) {
    // go through the PAF file
    igzstream sxs_in1(sxs_file.c_str());
    if (!sxs_in1.good()) assert("SXS is not good!");
    uint64_t lines = std::count(std::istreambuf_iterator<char>(sxs_in1),
                                std::istreambuf_iterator<char>(), '\n');
    // .sxs format is multiline
    uint64_t n_aln = 0;
    std::string line;
    sxs_in1.close();
    sxs_in1.open(sxs_file.c_str());
    //assert(sxs_in.good());
    while (std::getline(sxs_in1, line)) {
        //std::cerr << "a line " << line << std::endl;
        if (line[0] == 'A') ++n_aln;
    }
    igzstream sxs_in(sxs_file.c_str());
    //assert(sxs_in.good());
    //std::cerr << "I see " << n_aln << " alignments " << std::endl;
#pragma omp parallel for schedule(dynamic)
    for (size_t i = 0; i < n_aln; ++i) {
        //auto& line = lines[i];
        //std::string line;
        sxs_t sxs;
#pragma omp critical (sxs_in)
        sxs.load(sxs_in);
        assert(sxs.good());
        size_t query_idx = seqidx.rank_of_seq_named(sxs.query_sequence_name);
        size_t query_len = seqidx.nth_seq_length(query_idx);
        size_t target_idx = seqidx.rank_of_seq_named(sxs.target_sequence_name);
        size_t target_len = seqidx.nth_seq_length(target_idx);
        bool q_rev = sxs.b_rev();
        size_t q_all_pos = 1 + (q_rev ? seqidx.pos_in_all_seqs(query_idx, sxs.query_end, false) - 1
                                : seqidx.pos_in_all_seqs(query_idx, sxs.query_start, false));
        size_t t_all_pos = 1 + seqidx.pos_in_all_seqs(target_idx, sxs.target_start, false);
        pos_t q_pos = make_pos_t(q_all_pos, q_rev);
        pos_t t_pos = make_pos_t(t_all_pos, false);
        for (auto& c : sxs.cigar) {
            switch (c.op) {
            case 'M':
            {
                pos_t q_pos_match_start = q_pos;
                pos_t t_pos_match_start = t_pos;
                uint64_t match_len = 0;
                auto add_match = [&](void) {
                    if (match_len && match_len >= min_match_len) {
                        if (is_rev(q_pos)) {
                            aln_iitree_builder.add({offset(q_pos)+1, offset(q_pos_match_start)+1, make_pos_t(offset(t_pos)-1, true)});
                            aln_iitree_builder.add({offset(t_pos_match_start), offset(t_pos), make_pos_t(offset(q_pos_match_start), true)});
                        } else {
                            aln_iitree_builder.add({offset(q_pos_match_start), offset(q_pos), t_pos_match_start});
                            aln_iitree_builder.add({offset(t_pos_match_start), offset(t_pos), q_pos_match_start});
                        }
                    }
                };
                for (size_t i = 0; i < c.len; ++i) {
                    if (seqidx.at_pos(q_pos) == seqidx.at_pos(t_pos)
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

/*
void filter_alignments(mmmulti::map<pos_t, aln_pos_t>& aln_mm,
                       mmmulti::map<pos_t, pos_t>& aln_filt_mm,
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
