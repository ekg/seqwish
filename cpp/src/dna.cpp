#include "dna.hpp"

namespace seqwish {

static const char dna_complement[256] = {'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', // 8
                                         'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', // 16
                                         'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', // 24
                                         'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', // 32
                                         'N', 'N', 'N', '$', '#', 'N', 'N', 'N', // 40 GCSA stop/start characters
                                         'N', 'N', 'N', 'N', 'N', '-', 'N', 'N', // 48
                                         'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', // 56
                                         'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', // 64
                                         'N', 'T', 'V', 'G', 'H', 'N', 'N', 'C', // 72
                                         'D', 'N', 'N', 'M', 'N', 'K', 'N', 'N', // 80
                                         'N', 'Q', 'Y', 'W', 'A', 'A', 'B', 'S', // 88
                                         'N', 'R', 'N', 'N', 'N', 'N', 'N', 'N', // 96
                                         'N', 't', 'v', 'g', 'h', 'N', 'N', 'c', // 104
                                         'd', 'N', 'N', 'm', 'N', 'k', 'n', 'N', // 112
                                         'N', 'q', 'y', 'w', 'a', 'a', 'b', 's', // 120
                                         'N', 'r', 'N', 'N', 'N', 'N', 'N', 'N', // 128
                                         'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', // 136
                                         'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', // 144
                                         'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', // 152
                                         'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', // 160
                                         'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', // 168
                                         'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', // 176
                                         'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', // 184
                                         'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', // 192
                                         'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', // 200
                                         'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', // 208
                                         'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', // 216
                                         'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', // 224
                                         'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', // 232
                                         'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', // 240
                                         'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', // 248
                                         'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N'};// 256

char dna_reverse_complement(const char& c) {
    return dna_complement[c];
}

std::string dna_reverse_complement(const std::string& seq) {
    std::string rc;
    rc.assign(seq.rbegin(), seq.rend());
    for (auto& c : rc) {
        c = dna_complement[c];
    }
    return rc;
}
    
void dna_reverse_complement_in_place(std::string& seq) {
    size_t swap_size = seq.size() / 2;
    for (size_t i = 0, j = seq.size() - 1; i < swap_size; i++, j--) {
        char tmp = seq[i];
        seq[i] = dna_complement[seq[j]];
        seq[j] = dna_complement[tmp];
    }
    
    if (seq.size() % 2) {
        seq[swap_size] = dna_complement[seq[swap_size]];
    }
}

}
