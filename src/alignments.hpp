#pragma once

#include <iostream>
#include <string>
#include <vector>
#include <thread>
#include <mutex>
#include <atomic>
#include "paf.hpp"
#include "sxs.hpp"
#include "mmmultimap.hpp"
#include "mmiitree.hpp"
#include "seqindex.hpp"
#include "gzstream.h"
#include "pos.hpp"

namespace seqwish {

void paf_worker(
    igzstream& paf_in,
    std::atomic<bool>& paf_more,
    std::mutex& paf_in_mutex,
    mmmulti::iitree<uint64_t, pos_t>& aln_iitree,
    const seqindex_t& seqidx,
    const uint64_t& min_match_len);

void unpack_paf_alignments(
    const std::string& paf_file,
    mmmulti::iitree<uint64_t, pos_t>& aln_iitree,
    const seqindex_t& seqidx,
    const uint64_t& min_match_len,
    const uint64_t& num_threads);


/*
void filter_alignments(mmmulti::map<pos_t, aln_pos_t>& aln_mm,
                       mmmulti::map<pos_t, pos_t>& aln_filt_mm,
                       uint64_t aln_min_length,
                       uint64_t aln_keep_n_longest,
                       seqindex_t& seqidx);
*/

}
