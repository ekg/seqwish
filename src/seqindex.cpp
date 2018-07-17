#include "seqindex.hpp"

namespace seqwish {

// load a FASTA or FASTQ file into a file with a name index mapping name -> offset and indexed with a CSA
// provide queries over this index that let us extract particular positions and subsequences
void seqindex_t::set_base_filename(const std::string& filename) {
    basefilename = filename;
    seqfilename = basefilename + ".sqq";
    seqidxfile = basefilename + ".sqi";
    seqnamefile = basefilename + ".sqi.seqnames.tmp"; // used during construction
}

void seqindex_t::build_index(const std::string& filename) {
    set_base_filename(filename);
    // read the file
    igzstream in(filename.c_str());
    std::ofstream seqnames(seqnamefile.c_str());
    std::ofstream seqout(seqfilename.c_str());
    std::vector<uint64_t> seqname_offset;
    std::vector<uint64_t> seq_offset;
    bool input_is_fasta=false, input_is_fastq=false;
    // look at the first character to determine if it's fastq or fasta
    std::string line;
    std::getline(in, line);
    if (line[0] == '>') {
        input_is_fasta = true;
    } else if (line[0] == '@') {
        input_is_fastq = true;
    } else {
        std::cerr << "unknown file format given to seqindex_t" << std::endl;
        assert(false);
    }
    size_t seq_bytes_written = 0;
    size_t seq_names_bytes_written = 0;
    while (in.good()) {
        seqname_offset.push_back(seq_names_bytes_written);
        seq_offset.push_back(seq_bytes_written);
        line[0] = '>';
        line = line.substr(0, line.find(" "));
        seqnames << line;
        seq_names_bytes_written += line.size();
        std::string seq;
        // get the sequence
        if (input_is_fasta) {
            while (std::getline(in, line)) {
                if (line[0] == '>') {
                    // this is the header of the next sequence
                    break;
                } else {
                    seq.append(line);
                }
            }
        } else if (input_is_fastq) {
            std::getline(in, seq); // sequence
            std::getline(in, line); // delimiter
            std::getline(in, line); // quality
            std::getline(in, line);
        }
        seqout << seq;
        // record where the sequence starts
        seq_bytes_written += seq.size();
    }
    in.close();
    // add the last value so we can get sequence length for the last sequence and name
    seq_offset.push_back(seq_bytes_written);
    seqname_offset.push_back(seq_names_bytes_written);
    seqnames.close();
    seqout.close();
    // save the count of sequences
    seq_count = seqname_offset.size()-1;
    // mark the seq name starts vector, adding a terminating mark
    sdsl::bit_vector seq_name_starts(seqname_offset.back()+1);
    for (size_t i = 0; i < seqname_offset.size(); ++i) {
        seq_name_starts[seqname_offset[i]] = 1;
    }
    // build the name index
    construct(seq_name_csa, seqnamefile, 1);
    // destroy the file
    std::remove(seqnamefile.c_str());
    // build the rest of the index
    sdsl::util::assign(seq_name_cbv, sdsl::sd_vector<>(seq_name_starts));
    sdsl::util::assign(seq_name_cbv_rank, sdsl::sd_vector<>::rank_1_type(&seq_name_cbv));
    sdsl::util::assign(seq_name_cbv_select, sdsl::sd_vector<>::select_1_type(&seq_name_cbv));
    sdsl::util::assign(seq_offset_civ, sdsl::dac_vector<>(seq_offset));
    //std::cerr << seq_offset_civ << std::endl;
    // validate
    // look up each sequence by name
}

size_t seqindex_t::save(sdsl::structure_tree_node* s, std::string name) {
    assert(seq_name_csa.size() && seq_name_cbv.size() && seq_offset_civ.size());
    sdsl::structure_tree_node* child = sdsl::structure_tree::add_child(s, name, sdsl::util::class_name(*this));
    // open the sdsl index
    std::ofstream out(seqidxfile.c_str());
    size_t written = 0;
    out << "seqidx"; written += 9;
    uint32_t version_buffer = OUTPUT_VERSION;
    out.write((char*) &version_buffer, sizeof(version_buffer));
    written += sdsl::write_member(seq_count, out, child, "seq_count");
    written += seq_name_csa.serialize(out, child, "seq_name_csa");
    written += seq_name_cbv.serialize(out, child, "seq_name_cbv");
    written += seq_name_cbv_rank.serialize(out, child, "seq_name_cbv_rank");
    written += seq_name_cbv_select.serialize(out, child, "seq_name_cbv_select");
    written += seq_offset_civ.serialize(out, child, "seq_offset_civ");
    out.close();
    seqfile.open(seqfilename); // open the seq file
    return written;
}

void seqindex_t::load(const std::string& filename) {
    set_base_filename(filename);
    std::ifstream in(seqidxfile.c_str());
    std::string magic;
    in.read((char*)magic.c_str(), 6);
    uint32_t version;
    in.read((char*) &version, sizeof(version));
    assert(version == OUTPUT_VERSION);
    sdsl::read_member(seq_count, in);
    seq_name_csa.load(in);
    seq_name_cbv.load(in);
    seq_name_cbv_rank.load(in);
    seq_name_cbv_select.load(in);
    seq_offset_civ.load(in);
    in.close(); // close the sdsl index input
    seqfile.open(seqfilename); // open the seq file
}

void seqindex_t::to_fasta(std::ostream& out, size_t linewidth) {
    // extract the sequence names
    for (size_t i = 1; i < seq_count+1; ++i) {
        auto name = nth_name(i);
        out << ">" << name << std::endl;
        // pad sequence
        size_t seq_length = nth_seq_length(i);
        // for chunk of 80 up to sequence length
        //out << subseq(name, 0, linewidth) << std::endl;
        for (size_t j = 0; j < seq_length; j += linewidth) {
            out << subseq(name, j, std::min(linewidth, seq_length - j)) << std::endl;
        }
    }
    // iterate through the sequence names and extract the sequences
}

std::string seqindex_t::nth_name(size_t n) {
    // get the extents from our seq name dictionary
    size_t begin = seq_name_cbv_select(n)+1; // step past '>' delimiter
    size_t end = seq_name_cbv_select(n+1)-1;
    std::string name = sdsl::extract(seq_name_csa, begin, end);
    return name;
}

size_t seqindex_t::rank_of_seq_named(const std::string& name) {
    std::string query = ">" + name;
    //std::cerr << query << std::endl;
    auto occs = locate(seq_name_csa, query);
    //std::cerr << "occurs " << occs << std::endl;
    assert(occs.size() == 1);
    return seq_name_cbv_rank(occs[0])+1;
}

size_t seqindex_t::nth_seq_length(size_t n) {
    //std::cerr << "trying for "  << n << std::endl;
    return seq_offset_civ[n]-seq_offset_civ[n-1];
}

size_t seqindex_t::nth_seq_offset(size_t n) {
    return seq_offset_civ[n-1];
}

std::string seqindex_t::seq(const std::string& name) {
    return subseq(name, 0, nth_seq_length(rank_of_seq_named(name)));
}

std::string seqindex_t::subseq(const std::string& name, size_t pos, size_t count) {
    size_t n = rank_of_seq_named(name);
    return subseq(n, pos, count);
}

std::string seqindex_t::subseq(size_t n, size_t pos, size_t count) {
    return subseq(nth_seq_offset(n)+pos, count);
}

std::string seqindex_t::subseq(size_t pos, size_t count) {
    char s[count];
    seqfile.seekg(pos);
    seqfile.read(s, count);
    return std::string(s, count);
}

size_t seqindex_t::pos_in_all_seqs(const std::string& name, size_t pos, bool is_rev) {
    return pos_in_all_seqs(rank_of_seq_named(name), pos, is_rev);
}

size_t seqindex_t::pos_in_all_seqs(size_t n, size_t pos, bool is_rev) {
    return nth_seq_offset(n) + (is_rev ? nth_seq_length(n)-pos : pos);
}

size_t seqindex_t::seq_length(void) {
    return seq_offset_civ[seq_offset_civ.size()-1];
}

char seqindex_t::at(size_t pos) {
    char c;
    seqfile.seekg(pos);
    seqfile.read(&c, 1);
    return c;
}

char seqindex_t::at_pos(pos_t pos) {
    char c = at(offset(pos));
    if (is_rev(pos)) {
        c = dna_complement(c);
    }
    return c;
}

size_t seqindex_t::n_seqs(void) {
    return seq_count;
}

}
