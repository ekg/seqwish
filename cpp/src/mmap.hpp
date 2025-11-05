#ifndef MMAP_H_INCLUDED
#define MMAP_H_INCLUDED

#include <iostream>
#include <cstdio>
#include <string>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <fcntl.h>
#include <unistd.h>
#include <string>
#include <cassert>
#include "seqwish_rs.h"

namespace seqwish {

inline size_t mmap_open(const std::string& filename, char*& buf, int& fd) {
    char* buf_tmp = nullptr;
    int fd_tmp = -1;
    size_t size = ::mmap_open_rust(filename.c_str(), &buf_tmp, &fd_tmp);

    if (size == 0) {
        // Error occurred
        assert(false);
    }

    buf = buf_tmp;
    fd = fd_tmp;
    return size;
}

inline void mmap_close(char*& buf, int& fd, size_t fsize) {
    ::mmap_close_rust(buf, fd, fsize);
    buf = nullptr;
    fd = 0;
}

}

#endif
