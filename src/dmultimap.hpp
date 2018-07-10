#ifndef SEQWISH_HPP_INCLUDED
#define SEQWISH_HPP_INCLUDED

#include "bsort.hpp"
#include <iostream>
#include <string>
#include "sdsl/bit_vectors.hpp"

namespace seqwish {

template <typename K, typename V> class dmultimap {
    std::fstream file;
public:
    // constructor
    dmultimap(void) { }
    ~dmultimap(void) { }
    // close/open backing file
    bool open(const std::string& filename) {
        file.open(filename.c_str());
    }
    bool close(void) {
    }
    // push bacc (write to end of backing file)
    void push_back(const K& k, const V& v) {
    }

    // sort
    // index
    // load from base file name
};

}

#endif
