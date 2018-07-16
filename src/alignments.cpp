#include "alignments.hpp"

namespace seqwish {

void unpack_alignments(const std::string& paf_file,
                       dmultimap<int64_t, int64_t>& aln_mm,
                       const seqindex_t& seqidx) {
    // go through the PAF file
    std::ifstream paf_in(paf_file.c_str());
    std::string line;
    while (std::getline(paf_in, line)) {
        paf_row_t pafrow(line);
        // transform the alignments into a series of pairs b_i, b_j
        // that point into the seqindex
        //pafrow.cigar;
    }

}                      

}
