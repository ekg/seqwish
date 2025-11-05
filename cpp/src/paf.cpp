#include "paf.hpp"
#include "seqwish_rs.h"

// PAF row parsing now implemented in Rust (seqwish-rs/src/paf.rs)
// This file only provides utility functions

namespace seqwish {

std::ostream& operator<<(std::ostream& out, const paf_row_t& pafrow) {
    out << pafrow.query_sequence_name << "\t"
        << pafrow.query_sequence_length << "\t"
        << pafrow.query_start << "\t"
        << pafrow.query_end << "\t"
        << (pafrow.query_target_same_strand?"+":"-") << "\t"
        << pafrow.target_sequence_name << "\t"
        << pafrow.target_sequence_length << "\t"
        << pafrow.target_start << "\t"
        << pafrow.target_end << "\t"
        << pafrow.num_matches << "\t"
        << pafrow.alignment_block_length << "\t"
        << pafrow.mapping_quality << "\t"
        << "cg:Z:" << cigar_to_string(pafrow.cigar);
    return out;
}

void dump_paf_alignments(const std::string& filename) {
    std::ifstream in(filename.c_str());
    std::string line;
    while (std::getline(in, line)) {
        paf_row_t pafrow(line);
        std::cout << pafrow << std::endl;
    }
}

// Callback helper for parse_paf_spec FFI
extern "C" void paf_spec_callback(void* user_data, const char* filename, uint64_t weight) {
    auto* parsed = static_cast<std::vector<std::pair<std::string, uint64_t>>*>(user_data);
    parsed->push_back(std::make_pair(std::string(filename), weight));
}

std::vector<std::pair<std::string, uint64_t>> parse_paf_spec(const std::string& spec) {
    std::vector<std::pair<std::string, uint64_t>> parsed;
    ::parse_paf_spec(spec.c_str(), &parsed, paf_spec_callback);
    return parsed;
}

}
