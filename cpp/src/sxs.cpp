#include "sxs.hpp"

// SXS parsing now implemented in Rust (seqwish-rs/src/sxs.rs)
// This file only provides utility functions

namespace seqwish {

std::ostream& operator<<(std::ostream& out, const sxs_t& aln) {
    out << "A" << "\t" << aln.target_sequence_name << "\t" << aln.query_sequence_name << "\n"
        << "I" << "\t" << aln.target_start << "\t" << aln.target_end << "\t" << aln.query_start << "\t" << aln.query_end << "\n"
        << "M" << "\t" << aln.num_matches << "\n"
        << "C" << "\t" << cigar_to_string(aln.cigar) << "\n"
        << "Q" << "\t" << aln.mapping_quality;
    return out;
}

void dump_sxs_alignments(const std::string& filename) {
    std::ifstream in(filename.c_str());
    while (in.good()) {
        sxs_t aln(in);
        if (aln.good()) {
            std::cout << aln << std::endl;
        }
    }
}

}
