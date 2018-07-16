#include "paf.hpp"
#include "tokenize.hpp"

namespace seqwish {

PAFrow::PAFrow(const std::string& line) {
    std::vector<std::string> fields;
    tokenize(line, fields, " \t");
    query_sequence_name = fields[0];
    query_sequence_length = std::stol(fields[1]);
    query_start = std::stol(fields[2]);
    query_end = std::stol(fields[3]);
    query_target_same_strand = (fields[4] == "+");
    target_sequence_name = fields[5];
    target_sequence_length = std::stol(fields[6]);
    target_start = std::stol(fields[7]);
    target_end = std::stol(fields[8]);
    bases_in_mapping = std::stol(fields[9]);
    mapping_quality = std::stoi(fields[10]);
    // find the cigar in the last fields
    for (size_t i = 11; i < fields.size(); ++i) {
        // cg:Z:
        auto& f = fields[i];
        //std::string::size_type n;
        auto n = f.find("cg:Z:");
        if (n == 0) {
            cigar = cigar_from_string(f.substr(5));
            break;
        }
    }
}

std::ostream& operator<<(std::ostream& out, const PAFrow& pafrow) {
    out << pafrow.query_sequence_name << "\t"
        << pafrow.query_sequence_length << "\t"
        << pafrow.query_start << "\t"
        << pafrow.query_end << "\t"
        << pafrow.query_target_same_strand << "\t"
        << pafrow.target_sequence_name << "\t"
        << pafrow.target_sequence_length << "\t"
        << pafrow.target_start << "\t"
        << pafrow.target_end << "\t"
        << pafrow.bases_in_mapping << "\t"
        << pafrow.mapping_quality << "\t"
        << "cg:Z:" << cigar_to_string(pafrow.cigar);
    return out;
}

void parse_alignments(const std::string& filename) {
    std::ifstream in(filename.c_str());
    std::string line;
    while (std::getline(in, line)) {
        PAFrow pafrow(line);
        std::cout << pafrow << std::endl;
    }
}

}
