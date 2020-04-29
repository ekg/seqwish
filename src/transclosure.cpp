#include "transclosure.hpp"

namespace seqwish {

void extend_range(const uint64_t& s_pos,
                  const pos_t& q_pos,
                  std::map<pos_t, range_t>& range_buffer,
                  const seqindex_t& seqidx,
                  mmmulti::iitree<uint64_t, pos_t>& node_iitree,
                  mmmulti::iitree<uint64_t, pos_t>& path_iitree) {
    // find a position in the map that we can add onto
    // it must match position and orientation
    pos_t q_last_pos = q_pos;
    decr_pos(q_last_pos);
    auto f = range_buffer.find(q_last_pos);
    // if one doesn't exist, add the range
    if (f == range_buffer.end()) {
        range_buffer[q_pos] = {s_pos, s_pos+1};
    } else if ((!is_rev(q_pos) && seqidx.seq_start(offset(q_pos))) || (is_rev(q_pos) && seqidx.seq_start(offset(q_last_pos)))) {
        // flush the buffer we found, so we don't extend across node boundaries
        flush_range(f, node_iitree, path_iitree);
        range_buffer.erase(f);
        range_buffer[q_pos] = {s_pos, s_pos+1};
    } else {
        // if one does, check that it matches our extension,
        range_t x = f->second;
        if (x.end == s_pos) {
            // if so we expand its range and drop it back into the map at the new Q end pos
            range_buffer.erase(f);
            ++x.end; // increment the match length
            range_buffer[q_pos] = x; // and stash it
        } else {
            // if it doesn't, we store a new range
            range_buffer[q_pos] = {s_pos, s_pos+1};
        }
    }
}

void flush_ranges(const uint64_t& s_pos,
                  std::map<pos_t, range_t>& range_buffer,
                  mmmulti::iitree<uint64_t, pos_t>& node_iitree,
                  mmmulti::iitree<uint64_t, pos_t>& path_iitree) {
    // for each range, we're going to see if we've stepped more than one past the end
    // if we have, we'll write them out
    std::map<pos_t, range_t>::iterator it = range_buffer.begin();
    while (it != range_buffer.end()) {
        if (it->second.end != s_pos) {
            flush_range(it, node_iitree, path_iitree);
            it = range_buffer.erase(it);
        } else {
            ++it;
        }
    }
}

void flush_range(std::map<pos_t, range_t>::iterator it,
                 mmmulti::iitree<uint64_t, pos_t>& node_iitree,
                 mmmulti::iitree<uint64_t, pos_t>& path_iitree) {
    auto& range_in_s = it->second;
    uint64_t match_length, match_start_in_s, match_end_in_s, match_start_in_q, match_end_in_q;
    pos_t match_pos_in_q, match_pos_in_s, match_start_pos_in_q;
    pos_t match_end_pos_in_q = it->first;
    bool is_rev_match = is_rev(match_end_pos_in_q);
    if (!is_rev_match) {
        match_length = range_in_s.end - range_in_s.begin;
        match_start_in_s = range_in_s.begin;
        match_end_in_s = range_in_s.end;
        match_end_in_q = offset(match_end_pos_in_q) + 1;
        match_start_in_q = match_end_in_q - match_length;
        match_pos_in_s = make_pos_t(match_start_in_s, false);
        match_pos_in_q = make_pos_t(match_start_in_q, false);
    } else {
        match_length = range_in_s.end - range_in_s.begin;
        match_start_in_s = range_in_s.begin;
        match_end_in_s = range_in_s.end;
        match_end_in_q = offset(match_end_pos_in_q);
        decr_pos(match_end_pos_in_q, match_length);
        match_pos_in_s = make_pos_t(match_end_in_s-1, true);
        match_pos_in_q = make_pos_t(offset(match_end_pos_in_q)-1, true);
        match_start_in_q = match_end_in_q;
        match_end_in_q = offset(match_end_pos_in_q);
    }
    node_iitree.add(match_start_in_s, match_end_in_s, match_pos_in_q);
    path_iitree.add(match_start_in_q, match_end_in_q, match_pos_in_s);
}

// break the big range into its component ranges that we haven't already closed,
// breaking on sequence breaks
void for_each_fresh_range(const match_t& range,
                          const std::vector<bool>& seen_bv,
                          const std::function<void(match_t)>& lambda) {
    // walk range, breaking where we've seen it, emiting new ranges
    uint64_t p = range.start;
    pos_t t = range.pos;
    //std::cerr << "for_each_fresh_range " << range.start << "-" << range.end << " " << pos_to_string(range.pos) << std::endl;
    while (p < range.end) {
        // if we haven't seen p, start making a range
        //std::cerr << "looking at " << p << std::endl;
        if (seen_bv[p]) {
            ++p;
            incr_pos(t);
        } else {
            // otherwise, skip along
            uint64_t q = p;
            pos_t v = t;
            while (p < range.end && !seen_bv[p]) {
                ++p;
                incr_pos(t);
            }
            //std::cerr << "lambda\t" << q << " " << p << " " << pos_to_string(v) << std::endl;
            lambda({q, p, v});
        }
    }
}

void handle_range(match_t s,
                  atomicbitvector::atomic_bv_t& curr_bv,
                  overlap_atomic_queue_t& ovlp_q,
                  range_atomic_queue_t& todo_in) {
    bool all_set_there = true;
    pos_t n = s.pos;
    for (uint64_t i = s.start; i < s.end; ++i) {
        all_set_there = curr_bv.set(offset(n)) && all_set_there;
        incr_pos(n);
    }
    ovlp_q.push(std::make_pair(s, is_rev(s.pos)));
    if (!all_set_there) {
        auto item = std::make_pair(make_pos_t(offset(s.pos),is_rev(s.pos)), s.end - s.start);
        todo_in.push(item);
    }
}

void explore_overlaps(const match_t& b,
                      const std::vector<bool>& seen_bv,
                      atomicbitvector::atomic_bv_t& curr_bv,
                      mmmulti::iitree<uint64_t, pos_t>& aln_iitree,
                      overlap_atomic_queue_t& ovlp_q,
                      range_atomic_queue_t& todo_in) {
    //std::vector<size_t> o;
    aln_iitree.overlap(
        b.start, b.end,
        [&](const uint64_t& start,
            const uint64_t& end,
            const pos_t& pos) {
            match_t r = { start, end, pos };
            if (b.start > r.start) {
                uint64_t trim_from_start = b.start - r.start;
                r.start += trim_from_start;
                incr_pos(r.pos, trim_from_start);
            }
            if (r.end > b.end) {
                uint64_t trim_from_end = r.end - b.end;
                r.end -= trim_from_end;
            }
            assert(r.start < r.end);
            for_each_fresh_range(
                r,
                seen_bv,
                [&](match_t s) {
                    handle_range(s, curr_bv, ovlp_q, todo_in);
                });
        });
}

void write_graph_chunk(const seqindex_t& seqidx,
                       mmmulti::iitree<uint64_t, pos_t>& node_iitree,
                       mmmulti::iitree<uint64_t, pos_t>& path_iitree,
                       std::ofstream& seq_v_out,
                       std::map<pos_t, range_t>& range_buffer,
                       std::vector<std::pair<uint64_t, uint64_t>>* dsets_ptr,
                       uint64_t repeat_max,
                       uint64_t min_repeat_dist) {
    auto& dsets = *dsets_ptr;
    size_t seq_v_length = seq_v_out.tellp();
    uint64_t last_dset_id = std::numeric_limits<uint64_t>::max(); // ~inf
    char current_base = '\0';
    // determine if we've switched references
    // here we implement a count of the number of times we touch the current sequence
    std::map<uint64_t, uint64_t> seq_counts;
    std::map<uint64_t, pos_t> last_seq_pos;
    auto close_to_prev =
        [&last_seq_pos,
         &min_repeat_dist]
        (const uint64_t& seq_id,
         const pos_t& pos) {
            auto f = last_seq_pos.find(seq_id);
            if (f == last_seq_pos.end()) {
                return false;
            } else {
                if (min_repeat_dist > std::abs((int64_t)offset(pos) - (int64_t)offset(f->second))) {
                    return true;
                } else {
                    return false;
                }
            }
        };
    // run the closure for each dset, avoiding looping as configured
    std::map<uint64_t, std::vector<pos_t>> todos;
    std::string seq_out;
    auto flush_todos =
        [&](void) {
            for (auto& t : todos) {
                seq_out.push_back(current_base);
                ++seq_v_length;
                for (auto& pos : t.second) {
                    extend_range(seq_v_length-1, pos, range_buffer, seqidx, node_iitree, path_iitree);
                }
            }
        };
    for (auto& d : dsets) {
        const auto& curr_dset_id = d.first;
        const auto& curr_offset = d.second;
        char base = seqidx.at(curr_offset);
        // if we're on a new position
        if (curr_dset_id != last_dset_id) {
            if (repeat_max || min_repeat_dist) {
                // finish out todos stashed from repeat_max limitations
                flush_todos();
                todos.clear();
                // empty out our seq counts and last seq positions
                seq_counts.clear();
                last_seq_pos.clear();
            }
            // emit our new position
            current_base = base;
            seq_out.push_back(current_base);
            ++seq_v_length;
            flush_ranges(seq_v_length-1, range_buffer, node_iitree, path_iitree);
            last_dset_id = curr_dset_id;
        }
        pos_t curr_q_pos = make_pos_t(curr_offset, false);
        if (current_base != seqidx.at_pos(curr_q_pos)) {
            curr_q_pos = make_pos_t(curr_offset, true);
        }
        assert(current_base = seqidx.at_pos(curr_q_pos));
        uint64_t curr_seq_id = seqidx.seq_id_at(curr_offset);
        uint64_t curr_seq_count = 0;
        if ((min_repeat_dist != 0 && close_to_prev(curr_seq_id, curr_q_pos))
            || (repeat_max != 0 && seq_counts[curr_seq_id]+1 > repeat_max)) {
            curr_seq_count = ++seq_counts[curr_seq_id];
        } else if (repeat_max != 0 || min_repeat_dist != 0) {
            ++seq_counts[curr_seq_id];
        }
        if (curr_seq_count == 0) {
            extend_range(seq_v_length-1, curr_q_pos, range_buffer, seqidx, node_iitree, path_iitree);
        } else {
            todos[seq_counts[curr_seq_id]].push_back(curr_q_pos);
        }
        last_seq_pos[curr_seq_id] = curr_q_pos;
    }
    flush_todos(); // catch any todos we had hanging around
    seq_v_out << seq_out;
    delete dsets_ptr;
}


size_t compute_transitive_closures(
    const seqindex_t& seqidx,
    mmmulti::iitree<uint64_t, pos_t>& aln_iitree, // input alignment matches between query seqs
    const std::string& seq_v_file,
    mmmulti::iitree<uint64_t, pos_t>& node_iitree, // maps graph seq ranges to input seq ranges
    mmmulti::iitree<uint64_t, pos_t>& path_iitree, // maps input seq ranges to graph seq ranges
    uint64_t repeat_max,
    uint64_t min_repeat_dist,
    uint64_t transclose_batch_size, // size of a batch to collect for lock-free transitive closure
    bool show_progress,
    uint64_t num_threads,
    const std::chrono::time_point<std::chrono::steady_clock>& start_time) {
    // open the writers in the iitrees
    node_iitree.open_writer();
    path_iitree.open_writer();
    // open seq_v_file
    std::ofstream seq_v_out(seq_v_file.c_str());
    // remember the elements of Q we've seen
    //std::cerr << "seq_size " << seqidx.seq_length() << std::endl;
    std::vector<bool> q_seen_bv(seqidx.seq_length());
    //atomicbitvector::atomic_bv_t q_seen_bv(seqidx.seq_length());
    uint64_t input_seq_length = seqidx.seq_length();
    // a buffer of ranges to write into our iitree, arranged by range ending position in Q
    // we flush those intervals that don't get extended into the next position in S
    // this maps from a position in Q (our input seqs concatenated, offset and orientation)
    // to a range (start and length) in S (our graph sequence vector)
    // we are mapping from the /last/ position in the matched range, not the first
    std::map<pos_t, range_t> range_buffer;
    uint64_t bases_seen = 0;
    std::thread* graph_writer = nullptr;
    //uint64_t last_seq_id = seqidx.seq_id_at(0);
    // collect based on a seed chunk of a given length
    for (uint64_t i = 0; i < input_seq_length; ) {
        // scan our q_seen_bv to find our next start
        //std::cerr << "closing\t" << i << std::endl;
        while (i < input_seq_length && q_seen_bv[i]) ++i;
        //std::cerr << "scanned_to\t" << i << std::endl;
        if (i >= input_seq_length) break; // we're done!
        // where our chunk begins
        uint64_t chunk_start = i;
        // extend until we've got chunk_size unseen bases
        uint64_t bases_to_consider = 0;
        uint64_t chunk_end = chunk_start;
        while (bases_to_consider < transclose_batch_size && chunk_end < input_seq_length) {
            bases_to_consider += !q_seen_bv[chunk_end++];
        }
        // and where it ends (not past the end of the sequence)
        //chunk_end = std::min(input_seq_length, chunk_end); // chunk_start + transclose_batch_size);
        // collect ranges overlapping, per thread to avoid contention
        // bits of sequence we've seen during this union-find chunk
        atomicbitvector::atomic_bv_t q_curr_bv(seqidx.seq_length());
        // shared work queues for our threads
        auto todo_in_ptr = new range_atomic_queue_t;
        auto& todo_in = *todo_in_ptr;
        auto todo_out_ptr = new range_atomic_queue_t;
        auto& todo_out = *todo_out_ptr;
        auto ovlp_q_ptr = new overlap_atomic_queue_t;
        auto& ovlp_q = *ovlp_q_ptr;
        std::deque<std::pair<pos_t, uint64_t>> todo; // intermediate buffer for master thread
        std::vector<std::pair<match_t, bool>> ovlp; // written into by the master thread
#ifdef DEBUG_TRANSCLOSURE
        if (show_progress) std::cerr << "[seqwish::transclosure] " << std::fixed << std::showpoint << std::setprecision(3) << seconds_since(start_time) << " " << std::setprecision(2) << (double)bases_seen / (double)seqidx.seq_length() * 100 << "% " << chunk_start << "-" << chunk_end << " overlap_collect" << std::endl;
#endif
        // seed the initial ranges
        // the chunk range isn't an actual alignment, so we handle it differently
        for_each_fresh_range({chunk_start, chunk_end, 0}, q_seen_bv, [&](match_t b) {
                // the special case is handling ranges that have no matches
                // we need to close these even if they aren't matched to anything
                for (uint64_t j = b.start; j < b.end; ++j) {
                    assert(!q_seen_bv[j]);
                    q_curr_bv.set(j);
                }
                auto range = std::make_pair(make_pos_t(b.start, false), b.end - b.start);
                if (!todo_out.try_push(range)) {
                    todo.push_back(range);
                    //todo_seen.insert(range);
                }
            });
        std::atomic<bool> work_todo;
        std::vector<std::atomic<bool>> explorings(num_threads);
        work_todo.store(false);
        auto worker_lambda =
            [&](uint64_t tid) {
                //auto& ovlp = ovlps[tid];
                auto& exploring = explorings[tid];
                while (!work_todo.load()) {
                    std::this_thread::sleep_for(std::chrono::nanoseconds(1));
                }
                exploring.store(true);
                std::pair<pos_t, uint64_t> item;
                while (work_todo.load()) {
                    if (todo_out.try_pop(item)) {
                        exploring.store(true);
                        auto& pos = item.first;
                        auto& match_len = item.second;
                        uint64_t n = !is_rev(pos) ? offset(pos) : offset(pos) - match_len + 1;
                        uint64_t range_start = n;
                        uint64_t range_end = n + match_len;
                        explore_overlaps({range_start, range_end, pos},
                                         q_seen_bv,
                                         q_curr_bv,
                                         aln_iitree,
                                         ovlp_q,
                                         todo_in);
                    } else {
                        exploring.store(false);
                        std::this_thread::sleep_for(std::chrono::nanoseconds(1));
                    }
                }
                exploring.store(false);
            };
        // launch our threads to expand the overlap set in parallel
        std::vector<std::thread> workers; workers.reserve(num_threads);
        for (uint64_t t = 0; t < num_threads; ++t) {
            workers.emplace_back(worker_lambda, t);
        }
        // manage the threads
        uint64_t empty_iter_count = 0;
        auto still_exploring
            = [&explorings](void) {
                  bool ongoing = false;
                  for (auto& e : explorings) {
                      ongoing = ongoing || e.load();
                  }
                  return ongoing;
              };
        work_todo.store(true);
        while (!todo_in.was_empty() || !todo.empty() || !todo_out.was_empty() || !ovlp_q.was_empty() || still_exploring() || ++empty_iter_count < 1000) {
            std::this_thread::sleep_for(std::chrono::nanoseconds(1));
            // read from todo_in, into todo
            std::pair<pos_t, uint64_t> item;
            while (todo_in.try_pop(item)) {
                todo.push_back(item);
            }
            // then transfer to todo_out, until it's full
            while (!todo.empty()) {
                empty_iter_count = 0;
                item = todo.front();
                if (todo_out.try_push(item)) {
                    todo.pop_front();
                } else {
                    break;
                }
            }
            // collect our overlaps
            std::pair<match_t, bool> o;
            while (ovlp_q.try_pop(o)) {
                ovlp.push_back(o);
            }
        }
        //std::cerr << "telling threads to stop" << std::endl;
        work_todo.store(false);
        assert(todo.empty() && todo_in.was_empty() && todo_out.was_empty() && ovlp_q.was_empty() && !still_exploring());
        //std::cerr << "gonna join" << std::endl;
        for (uint64_t t = 0; t < num_threads; ++t) {
            workers[t].join();
        }
        delete todo_in_ptr;
        delete todo_out_ptr;
        delete ovlp_q_ptr;

#ifdef DEBUG_TRANSCLOSURE
        if (show_progress) std::cerr << "[seqwish::transclosure] " << std::fixed << std::showpoint << std::setprecision(3) << seconds_since(start_time) << " " << std::setprecision(2) << (double)bases_seen / (double)seqidx.seq_length() * 100 << "% " << chunk_start << "-" << chunk_end << " rank_build" << std::endl;
#endif
        // run the transclosure for this region using lock-free union find
        // convert the ranges into positions in the input sequence space
        uint64_t q_curr_bv_count = 0;
        //std::cerr << "q_subset_bv ";
        for (auto x : q_curr_bv) {
            ++q_curr_bv_count;
        }
        // use a rank support to make a dense mapping from the current bases to an integer range
        std::vector<uint64_t> q_curr_bv_vec; q_curr_bv_vec.reserve(q_curr_bv_count);
        for (auto p : q_curr_bv) {
            q_curr_bv_vec.push_back(p);
        }
        sdsl::bit_vector q_curr_bv_sdsl(seqidx.seq_length());
        for (auto p : q_curr_bv_vec) {
            q_curr_bv_sdsl[p] = 1;
        }
        sdsl::bit_vector::rank_1_type q_curr_rank;
        sdsl::util::assign(q_curr_rank, sdsl::bit_vector::rank_1_type(&q_curr_bv_sdsl));
        //q_curr_bv_vec.clear();
        // disjoint set structure
        std::vector<DisjointSets::Aint> q_sets_data(q_curr_bv_count);
        // this initializes everything
        auto disjoint_sets = DisjointSets(q_sets_data.data(), q_sets_data.size());
#ifdef DEBUG_TRANSCLOSURE
        if (show_progress) std::cerr << "[seqwish::transclosure] " << std::fixed << std::showpoint << std::setprecision(3) << seconds_since(start_time) << " " << std::setprecision(2) << (double)bases_seen / (double)seqidx.seq_length() * 100 << "% " << chunk_start << "-" << chunk_end << " parallel_union_find" << std::endl;
#endif
        paryfor::parallel_for<uint64_t>(
            0, ovlp.size(), num_threads, 10000,
            [&](uint64_t k) {
                auto& s = ovlp.at(k);
                auto& r = s.first;
                pos_t p = r.pos;
                for (uint64_t j = r.start; j != r.end; ++j) {
                    // unite both sides of the overlap
                    disjoint_sets.unite(q_curr_rank(j), q_curr_rank(offset(p)));
                    incr_pos(p);
                }
            });
        // now read out our transclosures
#ifdef DEBUG_TRANSCLOSURE
        if (show_progress) std::cerr << "[seqwish::transclosure] " << std::fixed << std::showpoint << std::setprecision(3) << seconds_since(start_time) << " " << std::setprecision(2) << (double)bases_seen / (double)seqidx.seq_length() * 100 << "% " << chunk_start << "-" << chunk_end << " dset_write" << std::endl;
#endif
        // maps from dset id to query base
        auto* dsets_ptr = new std::vector<std::pair<uint64_t, uint64_t>>(q_curr_bv_count);
        auto& dsets = *dsets_ptr;
        std::pair<uint64_t, uint64_t> max_pair = std::make_pair(std::numeric_limits<uint64_t>::max(), std::numeric_limits<uint64_t>::max());
        paryfor::parallel_for<uint64_t>(
            0, q_curr_bv_count, num_threads, 10000,
            [&](uint64_t j) {
                auto& p = q_curr_bv_vec[j];
                if (!q_seen_bv[p]) {
                    dsets[j] = std::make_pair(disjoint_sets.find(q_curr_rank(p)), p);
                } else {
                    dsets[j] = max_pair;
                }
            });
        //q_curr_bv_vec.clear();
        // remove excluded elements
        dsets.erase(std::remove_if(dsets.begin(), dsets.end(),
                                   [&max_pair](const std::pair<uint64_t, uint64_t>& x) {
                                       return x == max_pair;
                                   }),
                    dsets.end());
        // compress the dsets
#ifdef DEBUG_TRANSCLOSURE
        if (show_progress) std::cerr << "[seqwish::transclosure] " << std::fixed << std::showpoint << std::setprecision(3) << seconds_since(start_time) << " " << std::setprecision(2) << (double)bases_seen / (double)seqidx.seq_length() * 100 << "% " << chunk_start << "-" << chunk_end << " dset_compression" << std::endl;
#endif
        ips4o::parallel::sort(dsets.begin(), dsets.end());

        uint64_t c = 0;
        assert(dsets.size());
        uint64_t l = dsets.front().first;
        for (auto& d : dsets) {
            if (d.first != l) {
                ++c;
                l = d.first;
            }
            d.first = c;
        }
        /*
        for (auto& d : dsets) {
            std::cerr << "sdset\t" << d.first << "\t" << d.second << std::endl;
        }
        */
        // sort by the smallest starting position in each disjoint set
#ifdef DEBUG_TRANSCLOSURE
        if (show_progress) std::cerr << "[seqwish::transclosure] " << std::fixed << std::showpoint << std::setprecision(3) << seconds_since(start_time) << " " << std::setprecision(2) << (double)bases_seen / (double)seqidx.seq_length() * 100 << "% " << chunk_start << "-" << chunk_end << " dset_sort" << std::endl;
#endif
        std::vector<std::pair<uint64_t, uint64_t>> dsets_by_min_pos(c+1);
        for (uint64_t x = 0; x < c+1; ++x) {
            dsets_by_min_pos[x].second = x;
            dsets_by_min_pos[x].first = std::numeric_limits<uint64_t>::max();
        }
        for (auto& d : dsets) {
            uint64_t& minpos = dsets_by_min_pos[d.first].first;
            minpos = std::min(minpos, d.second);
        }
        ips4o::parallel::sort(dsets_by_min_pos.begin(), dsets_by_min_pos.end());
        /*
        for (auto& d : dsets_by_min_pos) {
            std::cerr << "sdset_min_pos\t" << d.second << "\t" << d.first << std::endl;
        }
        */
        // invert the naming
#ifdef DEBUG_TRANSCLOSURE
        if (show_progress) std::cerr << "[seqwish::transclosure] " << std::fixed << std::showpoint << std::setprecision(3) << seconds_since(start_time) << " " << std::setprecision(2) << (double)bases_seen / (double)seqidx.seq_length() * 100 << "% " << chunk_start << "-" << chunk_end << " dset_invert" << std::endl;
#endif
        std::vector<uint64_t> dset_names(c+1);
        uint64_t x = 0;
        for (auto& d : dsets_by_min_pos) {
            dset_names[d.second] = x++;
        }
        // rename sdsets and re-sort
        for (auto& d : dsets) {
            d.first = dset_names[d.first];
        }
        ips4o::parallel::sort(dsets.begin(), dsets.end());
        /*
        for (auto& d : dsets) {
            std::cerr << "sdset_rename\t" << d.first << "\t" << pos_to_string(d.second) << std::endl;
        }
        */
        // now, run the graph emission
#ifdef DEBUG_TRANSCLOSURE
        if (show_progress) std::cerr << "[seqwish::transclosure] " << std::fixed << std::showpoint << std::setprecision(3) << seconds_since(start_time) << " " << std::setprecision(2) << (double)bases_seen / (double)seqidx.seq_length() * 100 << "% " << chunk_start << "-" << chunk_end << " graph_emission" << std::endl;
#endif
        // mark q_seen_bv
        for (auto& d : dsets) {
            const auto& curr_offset = d.second;
            q_seen_bv[curr_offset] = 1;
            ++bases_seen;
        }
        // wait for completion of the last writer
        if (graph_writer != nullptr) {
            graph_writer->join();
            delete graph_writer;
        }
        // spawn the graph writer thread
        graph_writer = new std::thread(write_graph_chunk,
                                       std::ref(seqidx),
                                       std::ref(node_iitree),
                                       std::ref(path_iitree),
                                       std::ref(seq_v_out),
                                       std::ref(range_buffer),
                                       dsets_ptr,
                                       repeat_max,
                                       min_repeat_dist);
    }
    // clean up the last writer
    if (graph_writer != nullptr) {
        graph_writer->join();
        delete graph_writer;
    }
    // close the graph sequence vector
    size_t seq_bytes = seq_v_out.tellp();
    seq_v_out.close();
    flush_ranges(seq_bytes+1, range_buffer, node_iitree, path_iitree);
    assert(range_buffer.empty());
#ifdef DEBUG_TRANSCLOSURE
    if (show_progress) std::cerr << "[seqwish::transclosure] " << std::fixed << std::showpoint << std::setprecision(3) << seconds_since(start_time) << " " << std::setprecision(2) << (double)bases_seen / (double)seqidx.seq_length() * 100 << "% " << "building node_iitree and path_iitree indexes" << std::endl;
#endif
    // close writers
    node_iitree.close_writer();
    path_iitree.close_writer();
    // build node_mm and path_mm indexes
    node_iitree.index(num_threads);
    path_iitree.index(num_threads);
#ifdef DEBUG_TRANSCLOSURE
    if (show_progress) std::cerr << "[seqwish::transclosure] " << std::fixed << std::showpoint << std::setprecision(3) << seconds_since(start_time) << " " << std::setprecision(2) << (double)bases_seen / (double)seqidx.seq_length() * 100 << "% " << "done" << std::endl;
#endif
    return seq_bytes;
}

}
