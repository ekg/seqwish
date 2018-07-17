#include "transclosure.hpp"
#include "dmultimap.hpp"

namespace seqwish {

size_t compute_transitive_closures(
    seqindex_t& seqidx,
    dmultimap<uint64_t, pos_t>& aln_mm,
    const std::string& seq_v_file,
    dmultimap<uint64_t, pos_t>& node_mm,
    dmultimap<uint64_t, pos_t>& path_mm) {
    // open seq_v_file
    std::ofstream seq_v_out(seq_v_file.c_str());
    // remember the elements of Q we've seen
    sdsl::bit_vector q_seen_bv(seqidx.seq_length());
    //size_t seq_bytes = 0;
    // for each base in seqidx
    //   collect mapped bases
    //   mark our seen bitvector
    //   emit the base
    //   record entries in node_mm and path_mm
    uint64_t input_seq_length = seqidx.seq_length();
    for (uint64_t i = 1; i <= input_seq_length; ++i) {
        //std::cerr << "top " << i << std::endl;
        if (q_seen_bv[i-1]) continue;
        // write base
        seq_v_out << seqidx.at(i-1);
        size_t seq_v_length = seq_v_out.tellp();
        // mark current
        q_seen_bv[i-1] = 1;
        // emit current
        //node_mm.append(i, make_pos_t(i,false));
        //path_mm.append(i, make_pos_t(i,false));
        //std::cerr << "closure " << i << " to " << seq_v_length << std::endl;
        std::set<pos_t> todo;
        todo.insert(make_pos_t(i, false));
        while (!todo.empty()) {
            pos_t j = *todo.begin();
            /*
            std::cerr << seq_v_length << "\t"
                      << offset(j) << "\t"
                      << (is_rev(j)?"-":"+") << std::endl;
            */
            todo.erase(todo.begin());
            //std::cerr << todo.size() << std::endl;
            //assert(q_seen_bv[offset(j)-1]==1);
            node_mm.append(seq_v_length, j);
            path_mm.append(offset(j), make_pos_t(seq_v_length,is_rev(j)));
            aln_mm.for_values_of(offset(j), [&](const pos_t& pos) {
                    if (pos) {
                        uint64_t k = offset(pos);
                        //std::cerr << "looking " << k << std::endl;
                        if (k && !q_seen_bv[k-1]) {
                            //std::cerr << "closing " << k << std::endl;
                            q_seen_bv[k-1] = 1;
                            todo.insert(pos);
                        }
                    }
                });
            // get the values for each of the todo
            // and add them to todo if we haven't done'm
        }
    }
    // close the graph sequence vector
    size_t seq_bytes = seq_v_out.tellp();
    seq_v_out.close();
    // build node_mm and path_mm indexes
    node_mm.index(seq_bytes);
    path_mm.index(input_seq_length);
    return seq_bytes;
}

}
