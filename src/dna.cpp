#include "dna.hpp"

namespace seqwish {

char dna_complement(char c) {
    switch(c)
    {   
    case 'A':
        return 'T';
    case 'T':
        return 'A';
    case 'G':
        return 'C';
    case 'C':
        return 'G';
    default:
        return 'N';
    }   
}

}
