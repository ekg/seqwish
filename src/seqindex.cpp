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

void seqindex_t::build_index(const std::string& filename, const std::string& idxbasename) {
    set_base_filename(idxbasename);
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
        std::string seq_name = line.substr(0, line.find(" "));

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
        if (seq.empty()){
            std::cerr << "[seqwish] WARNING: input FASTA file contains empty sequences." << std::endl;
        } else {
            seqnames << seq_name << " ";
            seq_names_bytes_written += line.size() + 1;

            // force the sequence to be upper-case
            std::transform(seq.begin(), seq.end(), seq.begin(), [](char c) { return std::toupper(c); });
            seqout << seq;
            // record where the sequence starts
            seq_bytes_written += seq.size();
        }
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

    // check if there are duplicated sequence names
    std::ifstream seqnames_in(seqnamefile.c_str());
    bool duplicated_ids = false;
    while (std::getline(seqnames_in, line, ' ')) {
        //std::cout << line << " (" << locate(seq_name_csa, line).size() << ")" << std::endl;
        std::string query = line + " ";
        if(locate(seq_name_csa, query).size() > 1){
            duplicated_ids = true;
            break;
        }
    }
    seqnames_in.close();

    // destroy the file
    std::remove(seqnamefile.c_str());

    if (duplicated_ids){
        std::cerr << "[seqwish] ERROR: the input sequences have duplicated IDs." << std::endl;
        exit(1);
    }

    // build the rest of the index
    sdsl::util::assign(seq_name_cbv, sdsl::sd_vector<>(seq_name_starts));
    sdsl::util::assign(seq_name_cbv_rank, sdsl::sd_vector<>::rank_1_type(&seq_name_cbv));
    sdsl::util::assign(seq_name_cbv_select, sdsl::sd_vector<>::select_1_type(&seq_name_cbv));
    // mark the seq begin vector, adding a terminating mark
    sdsl::bit_vector seq_begin_bv(seq_offset.back()+1);
    for (size_t i = 0; i < seq_offset.size(); ++i) {
        seq_begin_bv[seq_offset[i]] = 1;
    }
    sdsl::util::assign(seq_begin_cbv, sdsl::sd_vector<>(seq_begin_bv));
    sdsl::util::assign(seq_begin_cbv_rank, sdsl::sd_vector<>::rank_1_type(&seq_begin_cbv));
    sdsl::util::assign(seq_begin_cbv_select, sdsl::sd_vector<>::select_1_type(&seq_begin_cbv));
    //std::cerr << seq_offset_civ << std::endl;
    // validate
    // look up each sequence by name
}

size_t seqindex_t::save(sdsl::structure_tree_node* s, std::string name) {
    //assert(seq_name_csa.size() && seq_name_cbv.size() && seq_offset_civ.size());
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
    written += seq_begin_cbv.serialize(out, child, "seq_begin_cbv");
    written += seq_begin_cbv_rank.serialize(out, child, "seq_begin_cbv_rank");
    written += seq_begin_cbv_select.serialize(out, child, "seq_begin_cbv_select");
    out.close();
    open_seq(seqfilename);
    return written;
}

void seqindex_t::remove_index_files(void) {
    std::remove(seqfilename.c_str());
    std::remove(seqidxfile.c_str());
}

void seqindex_t::open_seq(const std::string& filename) {
    if (seq_fd) return; //open
    assert(!filename.empty());
    // open in binary mode as we are reading from this interface
    seq_fd = open(filename.c_str(), O_RDWR);
    if (seq_fd == -1) {
        assert(false);
    }
    struct stat stats;
    if (-1 == fstat(seq_fd, &stats)) {
        assert(false);
    }
    seq_size = stats.st_size;
    if (!(seq_buf =
          (char*) mmap(NULL,
                       seq_size,
                       PROT_READ | PROT_WRITE,
                       MAP_SHARED,
                       seq_fd,
                       0))) {
        assert(false);
    }
    madvise((void*)seq_buf, seq_size, POSIX_MADV_WILLNEED | POSIX_MADV_SEQUENTIAL);
}

void seqindex_t::close_seq(void) {
    if (seq_buf) {
        munmap(seq_buf, seq_size);
        seq_buf = 0;
    }
    if (seq_fd) {
        close(seq_fd);
        seq_fd = 0;
    }
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
    seq_begin_cbv.load(in);
    seq_begin_cbv_rank.load(in);
    seq_begin_cbv_select.load(in);
    in.close(); // close the sdsl index input
    open_seq(filename);
}

void seqindex_t::to_fasta(std::ostream& out, size_t linewidth) const {
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

std::string seqindex_t::nth_name(size_t n) const {
    // get the extents from our seq name dictionary
    size_t begin = seq_name_cbv_select(n)+1; // step past '>' delimiter
    size_t end = seq_name_cbv_select(n+1)-2; // step back past added ' '
    std::string name = sdsl::extract(seq_name_csa, begin, end);
    return name;
}

size_t seqindex_t::rank_of_seq_named(const std::string& name) const {
    std::string query = ">" + name + " ";
    //std::cerr << query << std::endl;
    auto occs = locate(seq_name_csa, query);
    //std::cerr << "occurs " << occs << std::endl;
    assert(occs.size() == 1);
    return seq_name_cbv_rank(occs[0])+1;
}

size_t seqindex_t::nth_seq_length(size_t n) const {
    //std::cerr << "trying for "  << n << std::endl;
    return seq_begin_cbv_select(n+1)-seq_begin_cbv_select(n);
}

size_t seqindex_t::nth_seq_offset(size_t n) const {
    return seq_begin_cbv_select(n);
}

std::string seqindex_t::seq(const std::string& name) const {
    return subseq(name, 0, nth_seq_length(rank_of_seq_named(name)));
}

std::string seqindex_t::subseq(const std::string& name, size_t pos, size_t count) const {
    size_t n = rank_of_seq_named(name);
    return subseq(n, pos, count);
}

std::string seqindex_t::subseq(size_t n, size_t pos, size_t count) const {
    return subseq(nth_seq_offset(n)+pos, count);
}

std::string seqindex_t::subseq(size_t pos, size_t count) const {
    std::string s; s.resize(count);
    memcpy((void*)s.c_str(), &seq_buf[pos], count);
    return s;
}

size_t seqindex_t::pos_in_all_seqs(const std::string& name, size_t pos, bool is_rev) const {
    return pos_in_all_seqs(rank_of_seq_named(name), pos, is_rev);
}

size_t seqindex_t::pos_in_all_seqs(size_t n, size_t pos, bool is_rev) const {
    //std::cerr << "nth seq length " << nth_seq_length(n) << " offset " << nth_seq_offset(n) << std::endl;
    return nth_seq_offset(n) + (is_rev ? nth_seq_length(n)-1-pos : pos);
}

size_t seqindex_t::seq_length(void) const {
    return seq_begin_cbv.size()-1;
}

char seqindex_t::at(size_t pos) const {
    return seq_buf[pos];
}

char seqindex_t::at_pos(pos_t pos) const {
    // assumes 0-based pos
    char c = at(offset(pos));
    if (is_rev(pos)) {
        c = dna_reverse_complement(c);
    }
    return c;
}

size_t seqindex_t::n_seqs(void) const {
    return seq_count;
}

size_t seqindex_t::seq_id_at(size_t pos) const {
    return seq_begin_cbv_rank(pos+1);
}

bool seqindex_t::seq_start(size_t pos) const {
    return seq_begin_cbv[pos] == 1;
}

}
