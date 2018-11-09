#ifndef DMULTIMAP_HPP_INCLUDED
#define DMULTIMAP_HPP_INCLUDED

#include <iostream>
#include <string>
#include <sstream>
#include <functional>
#include <unordered_set>
#include "bsort.hpp"
#include "sdsl/bit_vectors.hpp"
#include "threads.hpp"

namespace seqwish {

/*
'dmultimap' is a disk-backed multimap where keys and values are stored
in a binary file. The key space is assumed to be numeric, but values
may be of arbitrary size.  To build the multimap we first append
key/value pairs.  To query the multimap we must first index it.  We
first sort by key using bsort.  Then we pad the key space so that we
have one entry per integer in the range [0, max(keys)], sorting again
to put the padding pairs in the right positions.  We record the key
space by marking a bitvector of length equal to max(keys) with 1 at
those positions corresponding to the first record of each key in the
sorted array.  We compress this bitvector and build select supports on
it We are now able to traverse the sorted array using select queries
on this bitvector.
*/

template <typename Key, typename Value> class dmultimap {

private:

    typedef struct { Key key; Value value; } Entry;
    std::ofstream writer;
    std::vector<std::ofstream> writers;
    char* reader;
    int reader_fd;
    std::string filename;
    std::string index_filename;
    bool sorted = false;
    // bsort parameters
    int char_start = 0;
    int char_stop = 255;
    int stack_size = 32;
    int cut_off = 32;
    size_t record_size = 0;
    // key information
    Key max_key = 0;
    // null key and value
    Key nullkey;
    Value nullvalue;
    // compressed bitvector marked at key starts
    sdsl::sd_vector<> key_cbv;
    //sdsl::sd_vector<>::rank_1_type key_cbv_rank;
    // select support for the key cbv
    sdsl::sd_vector<>::select_1_type key_cbv_select;
    bool indexed = false;
    uint32_t OUTPUT_VERSION = 1; // update as we change our format

    void init(void) {
        record_size = sizeof(Key) + sizeof(Value);
        nullkey = 0;
        for (size_t i = 0; i < sizeof(Value); ++i) {
            ((uint8_t*)&nullvalue)[i] = 0;
        }
    }

public:

    // constructor
    dmultimap(void) { init(); }

    dmultimap(const std::string& f) : filename(f) { init(); open_writers(f); }

    ~dmultimap(void) { close_writers(); }

    void set_base_filename(const std::string& f) {
        filename = f;
        index_filename = filename+".idx";
    }

    // load from base file name
    void load(const std::string& f) {
        open_reader();
        set_base_filename(f);
        std::ifstream in(index_filename.c_str());
        std::string magic;
        in.read((char*)magic.c_str(), 9);
        uint32_t version;
        in.read((char*) &version, sizeof(version));
        assert(version == OUTPUT_VERSION);
        size_t n_records, record_size_in_bytes;
        sdsl::read_member(record_size_in_bytes, in);
        assert(record_size_in_bytes == record_size);
        sdsl::read_member(n_records, in);
        assert(n_records == record_count());
        sdsl::read_member(max_key, in);
        assert(max_key == nth_key(n_records));
        key_cbv.load(in);
        key_cbv_select.load(in);
    }

    // save indexes
    size_t save(sdsl::structure_tree_node* s = NULL, std::string name = "") {
        assert(max_key && indexed);
        sdsl::structure_tree_node* child = sdsl::structure_tree::add_child(s, name, sdsl::util::class_name(*this));
        // open the sdsl index
        std::ofstream out(index_filename.c_str());
        size_t written = 0;
        out << "dmultimap"; written += 9;
        uint32_t version_buffer = OUTPUT_VERSION;
        out.write((char*) &version_buffer, sizeof(version_buffer));
        written += sdsl::write_member(record_size, out, child, "record_size");
        written += sdsl::write_member(record_count(), out, child, "record_count");
        written += sdsl::write_member(max_key, out, child, "max_key");
        written += key_cbv.serialize(out, child, "key_cbv");
        written += key_cbv_select.serialize(out, child, "key_cbv_select");
        out.close();
        return written;
    }

    // close/open backing file
    void open_main_writer(void) {
        if (writer.is_open()) {
            writer.seekp(0, std::ios_base::end); // seek to the end for appending
            return;
        }
        assert(!filename.empty());
        // open in binary append mode as that's how we write into the file
        writer.open(filename.c_str(), std::ios::binary | std::ios::app);
        if (writer.fail()) {
            throw std::ios_base::failure(std::strerror(errno));
        }
    }

    // per-thread writers
    void open_writers(const std::string& f) {
        set_base_filename(f);
        open_writers();
    }

    void open_writers(void) {
        assert(!filename.empty());
        writers.clear();
        writers.resize(get_thread_count());
        for (size_t i = 0; i < writers.size(); ++i) {
            auto& writer = writers[i];
            writer.open(writer_filename(i), std::ios::binary | std::ios::app);
            if (writer.fail()) {
                throw std::ios_base::failure(std::strerror(errno));
            }
        }
    }

    std::string writer_filename(size_t i) {
        std::stringstream wf;
        wf << filename << ".tmp_write" << "." << i;
        return wf.str();
    }

    void open_reader(void) {
        if (reader_fd) return; //open
        assert(!filename.empty());
        // open in binary mode as we are reading from this interface
        reader_fd = open(filename.c_str(), O_RDWR);
        if (reader_fd == -1) {
            assert(false);
        }
        struct stat stats;
        if (-1 == fstat(reader_fd, &stats)) {
            assert(false);
        }
        if (!(reader =
              (char*) mmap(NULL,
                           stats.st_size,
                           PROT_READ | PROT_WRITE,
                           MAP_SHARED,
                           reader_fd,
                           0))) {
            assert(false);
        }
        madvise((void*)reader, stats.st_size, POSIX_MADV_WILLNEED | POSIX_MADV_SEQUENTIAL);
    }

    std::ofstream& get_writer(void) {
        return writers[omp_get_thread_num()];
    }
    
    void sync_writers(void) {
        // close the temp writers and cat them onto the end of the main file
        open_main_writer();
        for (size_t i = 0; i < writers.size(); ++i) {
            writers[i].close();
            std::ifstream if_w(writer_filename(i), std::ios_base::binary);
            writer << if_w.rdbuf();
            if_w.close();
            std::remove(writer_filename(i).c_str());
        }
        writers.clear();
        writer.close();
    }

    void close_writers(void) {
        for (size_t i = 0; i < writers.size(); ++i) {
            std::remove(writer_filename(i).c_str());
        }
    }
    
    void close_reader(void) {
        if (reader) {
            size_t c = record_count();
            munmap(reader, c);
            reader = 0;
        }
        if (reader_fd) {
            close(reader_fd);
            reader_fd = 0;
        }
    }

    /// write the pair to end of backing file
    void append(const Key& k, const Value& v) {
        sorted = false; // assume we break the sort
        // write to the end of the file
        auto k_be = htobe64(k);
        auto& writer = get_writer();
        writer.write((char*)&k_be, sizeof(Key));
        writer.write((char*)&v, sizeof(Value));
    }

    /// get the record count
    size_t record_count(void) {
        int fd = open(filename.c_str(), O_RDWR);
        if (fd == -1) {
            assert(false);
        }
        struct stat stats;
        if (-1 == fstat(fd, &stats)) {
            assert(false);
        }
        assert(stats.st_size % record_size == 0); // must be even records
        size_t count = stats.st_size / record_size;
        return count;
    }

    /// sort the record in the backing file by key
    void sort(void) {
        if (sorted) return;
        //std::cerr << "sorting!" << std::endl;
        struct bsort::sort sort;
        if (-1==bsort::open_sort((char*)filename.c_str(), &sort)) {
            assert(false);
        }
        size_t key_size = sizeof(Key);
        bsort::radixify((unsigned char*)sort.buffer,
                        sort.size / record_size,
                        0,
                        char_start,
                        char_stop,
                        // size of the whole record
                        record_size,
                        // size of the key we sort by
                        record_size, // sort by the whole record so we can find unique values easily
                        stack_size,
                        cut_off);
        bsort::close_sort(&sort);
        sorted = true;
    }

    Entry read_entry(size_t i) {
        Entry e;
        //auto& reader = get_reader();
        memcpy(&e.key, &reader[i*record_size], sizeof(Key));
        memcpy(&e.value, &reader[i*record_size]+sizeof(Key), sizeof(Value));
        e.key = be64toh(e.key);
        return e;
    }

    // pad our key space with empty records so that we can query it directly with select operations
    void padsort(void) {
        close_reader();
        // blindly fill with a single key/value pair for each entity in the key space
        open_writers();
#pragma omp parallel for schedule(dynamic)
        for (size_t i = 1; i <= max_key; ++i) {
            append(i, nullvalue);
        }
        sync_writers();
        sort();
    }

    // index
    void index(Key new_max) {
        max_key = new_max;
        padsort();
        open_reader();
        size_t n_records = record_count();
        sdsl::bit_vector key_bv(n_records+1);
        // record the key starts
        Key last = nullkey, curr = nullkey;
        Value val = nullvalue;
        Entry entry;
        //reader.read((char*)&last, sizeof(Key));
        //key_bv[0] = 1;
        for (size_t i = 0; i < n_records; ++i) {
            entry = read_entry(i);
            curr = entry.key;
            val = entry.value;
            if (curr != last) {
                key_bv[i] = 1;
            }
            last = curr;
        }
        // the last key in the sort is our max key
        max_key = nth_key(n_records-1);
        key_bv[n_records] = 1; // sentinel
        // build the compressed bitvector
        sdsl::util::assign(key_cbv, sdsl::sd_vector<>(key_bv));
        key_bv.resize(0); // memory could be tight
        //sdsl::util::assign(key_cbv_rank, sdsl::sd_vector<>::rank_1_type(&key_cbv));
        // build the select supports on the key bitvector
        sdsl::util::assign(key_cbv_select, sdsl::sd_vector<>::select_1_type(&key_cbv));
        indexed = true;
    }

    Key nth_key(size_t n) {
        Entry e = read_entry(n);
        return e.key;
    }

    Value nth_value(size_t n) {
        Entry e = read_entry(n);
        return e.value;
    }

    void for_each_pair(const std::function<void(const Key&, const Value&)>& lambda) {
        open_reader(); // open or seek to beginning
        Key key;
        Value value;
        Entry entry;
        size_t n_records = record_count();
#pragma omp parallel for schedule(dynamic)
        for (size_t i = 0; i < n_records; ++i) {
            entry = read_entry(i);
            key = entry.key;
            value = entry.value;
            lambda(key, value);
        }
    }

    std::vector<Value> values(const Key& key) {
        std::vector<Value> values;
        for_values_of(key, [&values](const Value& v) { values.push_back(v); });
        return values;
    }

    std::vector<Value> unique_values(const Key& key) {
        std::vector<Value> values;
        for_unique_values_of(key, [&values](const Value& v) { values.push_back(v); });
        return values;
    }

    void for_unique_values_of(const Key& key, const std::function<void(const Value&)>& lambda) {
        // quirk: if we've sorted by the whole binary record,
        // then we can do a simple 'uniq' operation to get the unique values
        Value last = nullvalue;
        for_values_of(key, [this,&lambda,&last](const Value& value) {
                if (!is_null(value) && value != last) {
                    lambda(value);
                    last = value;
                }
            });
    }

    bool is_null(const Value& value) {
        for (size_t i = 0; i < sizeof(Value); ++i) {
            if (((uint8_t*)&value)[i] != 0) {
                return false;
            }
        }
        return true;
    }

    void for_values_of(const Key& key, const std::function<void(const Value&)>& lambda) {
        open_reader();
        size_t i = key_cbv_select(key);
        size_t j = key_cbv_select(key+1);
        for ( ; i < j; ++i) {
            Entry entry = read_entry(i);
            if (!is_null(entry.value)) {
                lambda(entry.value);
            }
        }
    }
};

}

#endif
