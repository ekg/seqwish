#ifndef SXS_HPP_INCLUDED
#define SXS_HPP_INCLUDED

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cassert>
#include "cigar.hpp"
#include "seqwish_rs.h"

namespace seqwish {

/// SXS alignment class - now backed by Rust implementation
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

    sxs_t(void) : query_sequence_length(0), query_start(0), query_end(0),
                  query_target_same_strand(true), target_sequence_length(0),
                  target_start(0), target_end(0), num_matches(0), mapping_quality(0) { }

    sxs_t(std::istream& in) {
        load(in);
    }

    void load(std::istream& in) {
        char c = in.get();
        if (in.eof()) return;
        // assert we have to start at the alignment
        assert(c == 'A');
        in.unget();

        // Collect lines for this alignment
        std::vector<std::string> lines;
        std::string line;
        std::getline(in, line);
        bool more = !line.empty();

        while (more) {
            lines.push_back(line);

            // check if we're to a new alignment
            if (!in.get(c)) {
                break;
            }
            in.unget();
            if (c == 'A') {
                more = false;
            } else {
                std::getline(in, line);
            }
        }

        // Parse using Rust
        std::vector<const char*> c_lines;
        for (const auto& l : lines) {
            c_lines.push_back(l.c_str());
        }

        struct SxsHandle* handle = ::sxs_parse_lines(c_lines.data(), c_lines.size());
        if (handle == nullptr) {
            // Failed to parse
            query_sequence_length = 0;
            return;
        }

        // Extract fields from Rust handle
        char* qname = ::sxs_query_sequence_name(handle);
        if (qname) {
            query_sequence_name = std::string(qname);
            ::temp_file_free_string(qname);
        }

        char* tname = ::sxs_target_sequence_name(handle);
        if (tname) {
            target_sequence_name = std::string(tname);
            ::temp_file_free_string(tname);
        }

        query_start = ::sxs_query_start(handle);
        query_end = ::sxs_query_end(handle);
        target_start = ::sxs_target_start(handle);
        target_end = ::sxs_target_end(handle);
        num_matches = ::sxs_num_matches(handle);
        mapping_quality = ::sxs_mapping_quality(handle);

        // Get CIGAR
        struct CigarHandle* cigar_handle = ::sxs_cigar(handle);
        if (cigar_handle) {
            cigar = cigar_from_handle(cigar_handle);
            ::cigar_free(cigar_handle);
        }

        // Free the handle
        ::sxs_free(handle);
    }

    friend std::ostream& operator<<(std::ostream& out, const sxs_t& aln);
};

void dump_sxs_alignments(const std::string& filename);

}

#endif
