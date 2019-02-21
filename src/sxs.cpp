#include "sxs.hpp"
#include "tokenize.hpp"

namespace seqwish {

sxs_t::sxs_t(std::istream& in) {
    load(in);
}

void sxs_t::load(std::istream& in) {
    char c = in.get();
    if (in.eof()) return;
    // assert we have to start at the alignment
    assert(c == 'A');
    in.unget();
    std::string line;
    std::getline(in, line);
    //std::cerr << "first line " << line << std::endl;
    bool more = !line.empty();
    while (more) {
        std::vector<std::string> fields;
        tokenize(line, fields, " \t");
        char c = fields[0][0];
        switch (c) {
        case 'A':
            target_sequence_name = fields[1];
            query_sequence_name = fields[2];
            break;
        case 'I':
            target_start = atoll(fields[1].c_str());
            target_end = atoll(fields[2].c_str());
            query_start = atoll(fields[3].c_str());
            query_end = atoll(fields[4].c_str());
            break;
        case 'M':
            num_matches = atoll(fields[1].c_str());
            break;
        case 'C':
            cigar = cigar_from_string(fields[1]);
            break;
        case 'T':
            // unhandled
            break;
        case 'Q':
            mapping_quality = std::stoi(fields[1]);
            break;
        default:
            break;
        }
        // check if we're to a new alignment
        if (!in.get(c)) {
            return;
        }
        in.unget();
        if (c == 'A') {
            more = false;
        } else {
            std::getline(in, line);
        }
    }
}

std::ostream& operator<<(std::ostream& out, const sxs_t& aln) {
    out << "A" << "\t" << aln.target_sequence_name << "\t" << aln.query_sequence_name << "\n"
        << "I" << "\t" << aln.target_start << "\t" << aln.target_end << "\t" << aln.query_start << "\t" << aln.query_end << "\n"
        << "M" << "\t" << aln.num_matches << "\n"
        << "C" << "\t" << cigar_to_string(aln.cigar) << "\n"
        << "Q" << "\t" << aln.mapping_quality;
    return out;
}

void dump_sxs_alignments(const std::string& filename) {
    std::ifstream in(filename.c_str());
    while (in.good()) {
        sxs_t aln(in);
        if (aln.good()) {
            std::cout << aln << std::endl;
        }
    }
}

}
