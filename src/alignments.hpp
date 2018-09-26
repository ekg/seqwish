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
                       dmultimap<pos_t, pos_t>& aln_mm,
                       seqindex_t& seqidx);

void filter_alignments(dmultimap<pos_t, aln_pos_t>& aln_mm,
                       dmultimap<pos_t, pos_t>& aln_filt_mm,
                       uint64_t aln_min_length,
                       uint64_t aln_keep_n_longest,
                       seqindex_t& seqidx);

}

#endif
