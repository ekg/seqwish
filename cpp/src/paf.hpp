#ifndef PAF_HPP_INCLUDED
#define PAF_HPP_INCLUDED

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cstdint>
#include "cigar.hpp"
#include "seqwish_rs.h"

namespace seqwish {

/// PAF row class - now backed by Rust implementation
class paf_row_t {
private:
    struct PafRowHandle* handle_;

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
    uint64_t alignment_block_length;
    uint16_t mapping_quality;
    cigar_t cigar;

    paf_row_t(const std::string& line) {
        // Parse using Rust
        handle_ = ::paf_row_parse(line.c_str());
        if (handle_ == nullptr) {
            // Failed to parse - leave fields at default values
            query_sequence_length = 0;
            return;
        }

        // Extract all fields from Rust handle
        char* qname = ::paf_row_query_sequence_name(handle_);
        if (qname) {
            query_sequence_name = std::string(qname);
            ::temp_file_free_string(qname);
        }

        char* tname = ::paf_row_target_sequence_name(handle_);
        if (tname) {
            target_sequence_name = std::string(tname);
            ::temp_file_free_string(tname);
        }

        query_sequence_length = ::paf_row_query_sequence_length(handle_);
        query_start = ::paf_row_query_start(handle_);
        query_end = ::paf_row_query_end(handle_);
        query_target_same_strand = ::paf_row_query_target_same_strand(handle_);
        target_sequence_length = ::paf_row_target_sequence_length(handle_);
        target_start = ::paf_row_target_start(handle_);
        target_end = ::paf_row_target_end(handle_);
        num_matches = ::paf_row_num_matches(handle_);
        alignment_block_length = ::paf_row_alignment_block_length(handle_);
        mapping_quality = ::paf_row_mapping_quality(handle_);

        // Get CIGAR
        struct CigarHandle* cigar_handle = ::paf_row_cigar(handle_);
        if (cigar_handle) {
            cigar = cigar_from_handle(cigar_handle);
            ::cigar_free(cigar_handle);
        }

        // Free the handle
        ::paf_row_free(handle_);
        handle_ = nullptr;
    }

    ~paf_row_t() {
        if (handle_ != nullptr) {
            ::paf_row_free(handle_);
        }
    }

    friend std::ostream& operator<<(std::ostream& out, const paf_row_t& pafrow);
};

void dump_paf_alignments(const std::string& filename);

std::vector<std::pair<std::string, uint64_t>> parse_paf_spec(const std::string& spec);

}

#endif
