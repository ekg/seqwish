#ifndef DNA_HPP_INCLUDED
#define DNA_HPP_INCLUDED

#include <string>

namespace seqwish {

char dna_reverse_complement(const char& c);
std::string dna_reverse_complement(const std::string& seq);
void dna_reverse_complement_in_place(std::string& seq);

}

#endif
