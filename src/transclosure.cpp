#include "transclosure.hpp"
#include "spinlock.hpp"
#include "dset64-gccAtomic.hpp"
#include "BooPHF.h"

namespace seqwish {

size_t compute_transitive_closures(
    seqindex_t& seqidx,
    range_pos_iitii& aln_iitree, // input alignment matches between query seqs
    const std::string& seq_v_file,
    mmmulti::iitree<uint64_t, pos_t>& node_iitree, // maps graph seq ranges to input seq ranges
    mmmulti::iitree<uint64_t, pos_t>& path_iitree, // maps input seq ranges to graph seq ranges
    uint64_t repeat_max,
    uint64_t min_transclose_len,
    uint64_t transclose_batch_size) { // size of a batch to collect for lock-free transitive closure
    // open seq_v_file
    std::ofstream seq_v_out(seq_v_file.c_str());
    // remember the elements of Q we've seen
    sdsl::bit_vector q_seen_bv(seqidx.seq_length());
    uint64_t input_seq_length = seqidx.seq_length();
    // a buffer of ranges to write into our iitree, arranged by range ending position in Q
    // we flush those intervals that don't get extended into the next position in S
    // this maps from a position in Q (our input seqs concatenated, offset and orientation)
    // to a range (start and length) in S (our graph sequence vector)
    // we are mapping from the /last/ position in the matched range, not the first
    std::map<pos_t, std::pair<uint64_t, uint64_t>> range_buffer;
    // here we try to find a growing range to extend
    auto extend_range = [&](const uint64_t& s_pos, const pos_t& q_pos) {
        // find a position in the map that we can add onto
        // it must match position and orientation
        pos_t q_last_pos = q_pos;
        decr_pos(q_last_pos);
        auto f = range_buffer.find(q_last_pos);
        // if one doesn't exist, add the range
        if (f == range_buffer.end()) {
            range_buffer[q_pos] = std::make_pair(s_pos, 1);
        } else {
            // if one does, check that it matches our extension,
            std::pair<uint64_t, uint64_t> x = f->second;
            if (x.first + x.second == s_pos) {
                // if so we expand its range and drop it back into the map at the new Q end pos
                range_buffer.erase(f);
                ++x.second; // increment the match length
                range_buffer[q_pos] = x; // and stash it
            } else {
                // if it doesn't, we store a new range
                range_buffer[q_pos] = std::make_pair(s_pos, 1);
            }
        }
    };
    auto flush_ranges = [&](const uint64_t& s_pos) {
        // for each range, we're going to see if we've stepped more than one past the end
        // if we have, we'll write them out
        std::map<pos_t, std::pair<uint64_t, uint64_t>>::iterator it = range_buffer.begin();
        while (it != range_buffer.end()) {
            if (it->second.first + it->second.second != s_pos) {
                // flush
                uint64_t match_length = it->second.second;
                uint64_t match_start_in_s = it->second.first;
                uint64_t match_end_in_s = match_start_in_s + match_length;
                pos_t match_end_pos_in_q = it->first;
                bool is_rev_match = is_rev(match_end_pos_in_q);
                // TODO appreciate why you're doing this here, presumably it's because of how we're running the transclosure
                if (!is_rev_match) incr_pos(match_end_pos_in_q, 1);
                //incr_pos(match_end_pos_in_q, 1);
                pos_t match_start_pos_in_q = match_end_pos_in_q;
                decr_pos(match_start_pos_in_q, match_length);
                uint64_t match_end_in_q = offset(match_end_pos_in_q);
                uint64_t match_start_in_q = offset(match_start_pos_in_q);
                pos_t match_pos_in_s = make_pos_t(match_start_in_s, is_rev_match);
                pos_t match_pos_in_q = make_pos_t(match_start_in_q, is_rev_match);
                if (is_rev_match) {
                    // go from transclosure model to the same pattern we have in the alignment iitree
                    // 1-based half open intervals, positions map to start and orientation in S and Q
                    std::swap(match_start_in_q, match_end_in_q);
                    incr_pos(match_pos_in_q, 1);
                    decr_pos(match_pos_in_s, match_length - 1);
                }
                node_iitree.add(match_start_in_s, match_end_in_s, match_pos_in_q);
                path_iitree.add(match_start_in_q, match_end_in_q, match_pos_in_s);
                it = range_buffer.erase(it);
            } else {
                ++it;
            }
        }
    };
    auto for_each_fresh_range = [&q_seen_bv](const range_pos_t& range,
                                             const std::function<void(range_pos_t)>& lambda) {
        // walk range, breaking where we've seen it, emiting new ranges
        uint64_t p = range.start;
        while (p < range.end) {
            // if we haven't seen p, start making a range
            if (q_seen_bv[p-1]) {
                ++p;
            } else {
                // otherwise, skip along
                uint64_t q = p;
                while (p < range.end && !q_seen_bv[p-1]) ++p;
                // recalculate the pos for our range
                pos_t t = range.pos;
                if (is_rev(t)) {
                    incr_pos(t, range.end - p);
                } else {
                    incr_pos(t, q - range.start);
                }
                lambda({q, p, t});
            }
        }
    };
    
    // range compressed model
    // accumulate pairs of ranges in S (output seq of graph) and Q (input seqs concatenated)
    // we flush when we stop extending
    // we determine when we stop extending when we have stepped a bp and broke our range extension
    //auto range_pos_hash = [](const range_pos_t& rp) { return std::hash<uint64_t>(rp.start) };
    uint64_t last_seq_id = seqidx.seq_id_at(1);
    for (uint64_t i = 1; i <= input_seq_length; i+=transclose_batch_size) {
        // collect ranges overlapping
        std::vector<std::pair<range_pos_t, bool>> ovlp;
        // complete our collection (todo: in parallel)
        std::set<std::tuple<uint64_t, uint64_t, pos_t>> seen;
        std::set<std::pair<pos_t, uint64_t>> todo;
        uint64_t chunk_start = i;
        uint64_t chunk_end = i + transclose_batch_size;
        for_each_fresh_range({chunk_start, chunk_end, 0}, [&](range_pos_t b) {
                for (auto& r : aln_iitree.overlap(b.start, b.end)) {
                    // XXX TODO break into ranges where we haven't already closed
                    for_each_fresh_range(r, [&](range_pos_t s) {
                            seen.insert({s.start, s.end, s.pos});
                            if (chunk_start > s.start) {
                                uint64_t trim_from_start = chunk_start - s.start;
                                s.start += trim_from_start;
                                incr_pos(s.pos, trim_from_start);
                            }
                            if (s.end > chunk_end) {
                                uint64_t trim_from_end = s.end - chunk_end;
                                s.end -= trim_from_end;
                            }
                            ovlp.push_back(std::make_pair(s, false));
                            todo.insert(std::make_pair(s.pos, s.end - s.start));
                        });
                }
            });
        while (!todo.empty()) {
            pos_t j = todo.begin()->first;
            uint64_t match_len = todo.begin()->second;
            todo.erase(todo.begin());
            // get the n and n+match_len on the forward strand
            uint64_t n = !is_rev(j) ? offset(j) : offset(j) - match_len;
            uint64_t range_start = n;
            uint64_t range_end = n + match_len;
            for (auto& r : aln_iitree.overlap(n, n+match_len)) {
                for_each_fresh_range(r, [&](range_pos_t s) {
                        // XXX TODO break into ranges where we haven't already closed
                        if (!seen.count({s.start, s.end, s.pos})) {
                            seen.insert({s.start, s.end, s.pos});
                            if (n > s.start) {
                                uint64_t trim_from_start = n - s.start;
                                s.start += trim_from_start;
                                incr_pos(s.pos, trim_from_start);
                            }
                            if (s.end > (n + match_len)) {
                                uint64_t trim_from_end = s.end - (n + match_len);
                                s.end -= trim_from_end;
                            }
                            // record the adjusted range and compute our flip relative to the tree of matches we're extending
                            ovlp.push_back(std::make_pair(s, is_rev(j)^is_rev(s.pos)));
                            todo.insert(std::make_pair(make_pos_t(offset(s.pos),is_rev(j)^is_rev(s.pos)), s.end - s.start));
                        }
                    });
            }
        }
        // print our overlaps
        std::cerr << "transc" << "\t" << i << "-" << i + transclose_batch_size << std::endl;
        for (auto& s : ovlp) {
            std::cerr << "ovlp" << "\t" << s.first.start << "-" << s.first.end << "\t" << offset(s.first.pos) << (is_rev(s.first.pos)?"-":"+") << "\t" << (s.second?"-":"+") << std::endl;
        }
        // run the transclosure for this region using lock-free union find
        
        // convert the ranges into positions in the input sequence space
        std::vector<pos_t> q_subset;
        for (auto& s : ovlp) {
            auto& r = s.first;
            pos_t p = r.pos;
            //bool flip = s.second;
            for (uint64_t j = r.start; j != r.end; ++j) {
                // XXX TODO skip if we've already seen it
                q_subset.push_back(make_pos_t(j, false));
                q_subset.push_back(make_pos_t(j, true));
                q_subset.push_back(p);
                q_subset.push_back(rev_pos_t(p));
                incr_pos(p);
            }
        }
        std::sort(q_subset.begin(), q_subset.end());
        q_subset.erase(std::unique(q_subset.begin(), q_subset.end()), q_subset.end());
        uint nthreads = get_thread_count();
        double gammaFactor = 2.0; // lowest bit/elem is achieved with gamma=1, higher values lead to larger mphf but faster construction/query
                                  // gamma = 2 is a good tradeoff (leads to approx 3.7 bits/key )
        // how to query: bphf.lookup(key) maps ids into [0..N)
        // make a dense mapping
        auto bphf = boomphf::mphf<pos_t, pos_t_hasher>(q_subset.size(),q_subset,nthreads,gammaFactor,false,false);
        // disjoint set structure
        std::vector<DisjointSets::Aint> q_sets_data(q_subset.size());
        // this initializes everything
        auto disjoint_sets = DisjointSets(q_sets_data.data(), q_sets_data.size());
#pragma omp parallel for
        for (uint64_t x = 0; x < q_subset.size(); ++x) {
            uint64_t j = offset(q_subset.at(x));
            // unite the forward and reverse strands for the given base
            disjoint_sets.unite(bphf.lookup(make_pos_t(j, false)), bphf.lookup(make_pos_t(j, true)));
        }
        // join our strands
#pragma omp parallel for
        for (uint64_t k = 0; k < ovlp.size(); ++k) {
            auto& s = ovlp.at(k);
            auto& r = s.first;
            pos_t p = r.pos;
            //bool flip = s.second;
            for (uint64_t j = r.start; j != r.end; ++j) {
                // XXX todo skip if we've already closed this base
                // unite both sides of the overlap
                disjoint_sets.unite(bphf.lookup(make_pos_t(j, false)), bphf.lookup(p));
                incr_pos(p);
            }
        }
        // now read out our transclosures
        std::vector<std::pair<uint64_t, pos_t>> dsets;
        for (uint64_t x = 0; x < q_subset.size(); ++x) {
            //uint64_t j = offset(q_subset.at(x));
            auto& p = q_subset.at(x);
            //disjoint_sets.unite(bphf.lookup(make_pos_t(j, false)), bphf.lookup(make_pos_t(j, true)));
            std::cerr << "dset\t" << pos_to_string(p) << "\t"
                      << disjoint_sets.find(bphf.lookup(p)) << std::endl;
            dsets.push_back(std::make_pair(disjoint_sets.find(bphf.lookup(p)), p));
        }
        std::sort(dsets.begin(), dsets.end());
        for (auto& d : dsets) {
            // run over the
            std::cerr << "sdset\t" << d.first << "\t" << pos_to_string(d.second) << std::endl;
        }

        // XXX TODO set our seen bv for all the things we're closing

//#pragma omp parallel for
        
        
        // we'll use this for random access to the probably discontiguous membership array
        

        // (todo: use this to initially verify that we close everything)
        
        // insert these integers into our disjoint set union find structure
        // (hopefully we can avoid recording the mappping)
        
        // in parallel, run over each ovlp, unioning the values in the disjoint set
        
        // sort by minimum position in each union of matched query positions
        
        // iterate over the result with our range compressor
        // writing out the bases in the same order that we would have with the deterministic transclosure
    }
    exit(1);
    // close the graph sequence vector
    size_t seq_bytes = seq_v_out.tellp();
    seq_v_out.close();
    flush_ranges(seq_bytes+2);
    assert(range_buffer.empty());
    // build node_mm and path_mm indexes
    node_iitree.index();
    path_iitree.index();
    return seq_bytes;
}

}
