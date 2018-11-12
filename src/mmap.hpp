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

namespace seqwish {

size_t mmap_open(const std::string& filename, char*& buf, int& fd);
void mmap_close(char*& buf, int& fd, size_t fsize);

}

#endif
