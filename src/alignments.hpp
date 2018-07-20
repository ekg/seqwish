#ifndef ALIGNMENTS_HPP_INCLUDED
#define ALIGNMENTS_HPP_INCLUDED

#include <iostream>
#include <string>
#include <vector>
#include "paf.hpp"
#include "dmultimap.hpp"
#include "seqindex.hpp"
#include "pos.hpp"

namespace seqwish {

void unpack_alignments(const std::string& paf_file,
                       seqwish::dmultimap<pos_t, aln_pos_t>& aln_mm,
                       seqindex_t& seqidx);

}

#endif
