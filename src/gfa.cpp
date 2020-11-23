#include "gfa.hpp"

namespace seqwish {

using namespace std::chrono_literals;

void emit_gfa(std::ostream& out,
              size_t graph_length,
              const std::string& seq_v_file,
              mmmulti::iitree<uint64_t, pos_t>& node_iitree,
              mmmulti::iitree<uint64_t, pos_t>& path_iitree,
              const sdsl::sd_vector<>& seq_id_cbv,
              const sdsl::sd_vector<>::rank_1_type& seq_id_cbv_rank,
              const sdsl::sd_vector<>::select_1_type& seq_id_cbv_select,
              seqindex_t& seqidx,
              mmmulti::set<std::pair<pos_t, pos_t>>& link_mmset,
              const uint64_t& num_threads) {

    out << "H" << "\t" << "VN:Z:1.0" << std::endl;
    int seq_v_fd = -1;
    char* seq_v_buf = nullptr;
    size_t seq_v_filesize = mmap_open(seq_v_file, seq_v_buf, seq_v_fd);

    auto show_links = [&](const pos_t& p) { std::cerr << pos_to_string(p) << " " << pos_to_string(make_pos_t(seq_id_cbv_rank(offset(p)), is_rev(p))) << ", "; };

    // write the nodes
    // these are delimited in the seq_v_file by the markers in seq_id_civ

    uint64_t n_nodes = seq_id_cbv_rank(seq_id_cbv.size()-1);

    // producer/consumer queues
    auto seq_todo_q_ptr = new atomic_queue::AtomicQueue2<uint64_t, 2 << 16>;
    auto& seq_todo_q = *seq_todo_q_ptr;
    auto seq_done_q_ptr = new atomic_queue::AtomicQueue2<std::pair<uint64_t, std::string*>, 2 << 16>;
    auto& seq_done_q = *seq_done_q_ptr;
    std::atomic<bool> work_todo;
    std::map<uint64_t, std::string*> node_records;

    auto worker_lambda =
        [&](void) {
            uint64_t id = 0;
            while (work_todo.load()) {
                if (seq_todo_q.try_pop(id)) {
                    //std::stringstream s;
                    size_t node_start = seq_id_cbv_select(id);
                    //size_t node_length = (id==n_nodes ? seq_id_cbv.size() : seq_id_cbv_select(id+1)) - node_start;
                    size_t node_length = seq_id_cbv_select(id+1) - node_start;
                    //std::cerr << id << " "  << node_start << " " << node_length << std::endl;
                    std::string* seq = new std::string;
                    seq->resize(node_length);
                    memcpy((void*)seq->c_str(), &seq_v_buf[node_start], node_length);
                    seq_done_q.push(std::make_pair(id, seq));
                } else {
                    std::this_thread::sleep_for(0.00001ns);
                }
            }
        };

    std::vector<std::thread> workers; workers.reserve(num_threads);
    work_todo.store(true);
    for (uint64_t t = 0; t < num_threads; ++t) {
        workers.emplace_back(worker_lambda);
    }
    uint64_t todo_id = 1;
    uint64_t done_id = 0;
    while (done_id < n_nodes) {
        // put ids in todo
        while (todo_id <= n_nodes && seq_todo_q.try_push(todo_id)) {
            ++todo_id;
        }
        // read from done queue
        std::pair<uint64_t, std::string*> item;
        while (seq_done_q.try_pop(item)) {
            node_records[item.first] = item.second;
        }
        if (node_records.size()) {
            auto b = node_records.begin();
            if (b->first == done_id+1) {
                //out << node_records.begin()->second << std::endl;
                out << "S" << "\t" << b->first << "\t" << *b->second << "\n";
                ++done_id;
                delete b->second;
                node_records.erase(b);
            }
        }
    }
    work_todo.store(false);
    for (uint64_t t = 0; t < num_threads; ++t) {
        workers[t].join();
    }

    auto print_link = [&out](const std::pair<pos_t, pos_t>& p) {
        auto& from = p.first;
        auto& to = p.second;
        if (from && to) {
            out << "L" << "\t"
                << offset(from) << "\t" << (is_rev(from)?"-":"+") << "\t"
                << offset(to) << "\t" << (is_rev(to)?"-":"+") << "\t"
                << "0M" << "\n";
        }
    };
    link_mmset.for_each_unique_value(print_link);

    // write the paths
    size_t num_seqs = seqidx.n_seqs();
    // should be made to run in parallel
    for (size_t i = 1; i <= num_seqs; ++i) {
        size_t j = seqidx.nth_seq_offset(i);
        size_t seq_len = seqidx.nth_seq_length(i);
        size_t k = j + seq_len;
        //std::cerr << seqidx.nth_name(i) << " " << seqidx.nth_seq_length(i) << " " << j << " " << k << std::endl;
        std::vector<pos_t> path_v;
        uint64_t seen_bp = 0;
        while (j < k) {
            uint64_t overlap_count = 0;
            uint64_t ovlp_start_in_q;
            uint64_t ovlp_end_in_q;
            pos_t pos_start_in_s;
            path_iitree.overlap(
                j, j+1,
                [&](const uint64_t& start,
                    const uint64_t& end,
                    const pos_t& pos) {
                    ++overlap_count;
                    ovlp_start_in_q = start;
                    ovlp_end_in_q = end;
                    pos_start_in_s = pos;
                });
            // each input base should only map one place in the graph
            if (overlap_count != 1) {
                std::cerr << "[seqwish::gfa] error: found " << overlap_count << " overlaps for seq " << seqidx.nth_name(i) << " idx " << i << " at j=" << j << " of " << k << std::endl;
                path_iitree.overlap(
                    j, j+1,
                    [&](const uint64_t& start,
                    const uint64_t& end,
                    const pos_t& pos) {
                        std::cerr << "ovlp_start_in_q = " << start << " "
                                  << "ovlp_end_in_q = " << end << " "
                                  << "pos_start_in_s = " << pos_to_string(pos) << std::endl;
                    });
                assert(false);
                exit(1);
            }
            bool match_is_rev = is_rev(pos_start_in_s);
            // iterate through the nodes in this range
            uint64_t length = ovlp_end_in_q - ovlp_start_in_q;
            //std::cerr << "overlap in path_iitree " << ovlp_start_in_q << ".." << ovlp_end_in_q << " " << pos_to_string(pos_start_in_s) << std::endl;
            pos_t q = make_pos_t(ovlp_start_in_q, false);
            pos_t p = pos_start_in_s;
            // validate the path
            // for each base in the range pointed to by the match, check that the input sequence we're processing matches the graph
            for (uint64_t k = 0; k < length; ++k) {
                if (seq_id_cbv[offset(p)]) {
                    uint64_t node_id = seq_id_cbv_rank(offset(p)+1);
                    //std::cerr << "got to node " << node_id << (match_is_rev ? "-" : "+") << std::endl;
                    path_v.push_back(make_pos_t(node_id, match_is_rev));
                }
                /*else if (match_is_rev && seq_id_cbv[offset(p)]) {
                    uint64_t node_id = seq_id_cbv_rank(offset(p)-1);
                    std::cerr << "got to node " << node_id << "-" << std::endl;
                    path_v.push_back(make_pos_t(node_id, match_is_rev));
                }
                */
                char c = seq_v_buf[offset(p)];
                if (is_rev(p)) c = dna_reverse_complement(c);
                //std::cerr << pos_to_string(q) << " -> " << pos_to_string(p) << " " << seqidx.at_pos(q) << " vs " << c << " " << seq_id_cbv[offset(p)] << std::endl;
                if (seqidx.at_pos(q) != c) {
                    std::cerr << "GRAPH BROKEN @ "
                        << seqidx.nth_name(i) << " " << pos_to_string(q) << " -> "
                        << pos_to_string(q) << std::endl;
                    assert(false);
                    exit(1); // for release builds
                }
                incr_pos(p, 1);
                incr_pos(q, 1);
            }
            seen_bp += length;
            j = ovlp_end_in_q;
        }
        if (seen_bp != seq_len) {
            std::cerr << "length for " << seqidx.nth_name(i) << ", expected " << seqidx.nth_seq_length(i) << " but got " << seen_bp << std::endl;
            assert(false);
            exit(1); // for release builds
        }
        std::stringstream cigarss;
        std::stringstream pathss;

        for (auto& p : path_v) {
            pathss << pos_to_string(p) << ",";
        }
        pathss.seekp(-1, pathss.cur); // trim the last ","
        cigarss << "*";
        pathss << '\t';
        //cigarss.seekp(-1, cigarss.cur);
        cigarss << "\n";
//#pragma omp critical (out)
        out << "P" << "\t"
            << seqidx.nth_name(i) << "\t"
            << pathss.str()
            << cigarss.str();
    }

    mmap_close(seq_v_buf, seq_v_fd, seq_v_filesize);

}

}
