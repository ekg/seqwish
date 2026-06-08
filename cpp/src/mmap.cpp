#include "mmap.hpp"

namespace seqwish {

size_t mmap_open(const std::string& filename, char*& buf, int& fd) {
    fd = -1;
    assert(!filename.empty());
    // open in binary mode as we are reading from this interface
    fd = open(filename.c_str(), O_RDWR);
    if (fd == -1) {
        assert(false);
    }
    struct stat stats;
    if (-1 == fstat(fd, &stats)) {
        assert(false);
    }
    size_t fsize = stats.st_size;
    if (!(buf =
          (char*) mmap(NULL,
                       fsize,
                       PROT_READ | PROT_WRITE,
                       MAP_SHARED,
                       fd,
                       0))) {
        assert(false);
    }
    madvise((void*)buf, fsize, POSIX_MADV_WILLNEED | POSIX_MADV_SEQUENTIAL);
    return fsize;
}

void mmap_close(char*& buf, int& fd, size_t fsize) {
    if (buf) {
        munmap(buf, fsize);
        buf = 0;
    }
    if (fd) {
        close(fd);
        fd = 0;
    }
}

}
