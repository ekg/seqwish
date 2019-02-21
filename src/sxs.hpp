#ifndef SXS_HPP_INCLUDED
#define SXS_HPP_INCLUDED

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cassert>
#include "cigar.hpp"

namespace seqwish {

class sxs_t {
public:
    std::string query_sequence_name;
    uint64_t query_sequence_length;
    uint64_t query_start;
    uint64_t query_end;
    bool query_target_same_strand;
    std::string target_sequence_name;
    uint64_t target_sequence_length;
    uint64_t target_start;
    uint64_t target_end;
    uint64_t num_matches;
    uint16_t mapping_quality;
    cigar_t cigar;
    bool good(void) { return query_sequence_name.size() > 0; }
    bool b_rev(void) { return query_start > query_end; }
    sxs_t(void) { }
    sxs_t(std::istream& in);
    void load(std::istream& in);
    friend std::ostream& operator<<(std::ostream& out, const sxs_t& aln);
};

void dump_sxs_alignments(const std::string& filename);

}

#endif
