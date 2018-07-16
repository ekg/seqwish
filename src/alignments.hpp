#ifndef ALIGNMENTS_HPP_INCLUDED
#define ALIGNMENTS_HPP_INCLUDED

#include <iostream>
#include <string>
#include <vector>
#include "paf.hpp"
#include "dmultimap.hpp"
#include "seqindex.hpp"

namespace seqwish {

void unpack_alignments(const std::string& paf_file,
                       seqwish::dmultimap<int64_t, int64_t>& aln_mm,
                       const seqindex_t& seqidx);

}

#endif
