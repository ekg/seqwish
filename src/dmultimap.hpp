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

    std::ofstream writer;
    std::vector<std::ofstream> writers;
    std::vector<std::ifstream> readers;
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
        open_readers();
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

    void open_readers(void) {
        if (readers.size()) {
            for (auto& reader : readers) {
                reader.seekg(0); // reset to start
            }
            return;
        } else {
            readers.resize(get_thread_count());
        }
        assert(!filename.empty());
        // open in binary mode as we are reading from this interface
        for (auto& reader : readers) {
            reader.open(filename, std::ios::binary);
            if (reader.fail()) {
                throw std::ios_base::failure(std::strerror(errno));
            }
        }
    }

    std::ifstream& get_reader(void) {
        return readers[omp_get_thread_num()];
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
    
    void close_readers(void) {
        readers.clear();
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
        open_readers();
        auto& reader = get_reader();
        auto pos = reader.tellg();
        reader.seekg(0, std::ios_base::end); // seek to the end
        assert(reader.tellg() % record_size == 0); // must be even records
        size_t count = reader.tellg() / record_size;
        reader.seekg(pos);
        return count;
    }

    /// sort the record in the backing file by key
    void sort(void) {
        if (sorted) return;
        //std::cerr << "sorting!" << std::endl;
        close_readers();
        sync_writers();
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
                        record_size,
                        key_size,
                        stack_size,
                        cut_off);
        bsort::close_sort(&sort);
        sorted = true;
    }

    Key read_key(void) {
        Key k;
        auto& reader = get_reader();
        reader.read((char*)&k, sizeof(Key));
        return be64toh(k); // sorting is big-endian
    }

    Value read_value(void) {
        Value v;
        auto& reader = get_reader();
        reader.read((char*)&v, sizeof(Value));
        return v;
    }

    // pad our key space so that we can query it directly with select operations
    void pad(void) {
        assert(sorted);
        open_readers();
        // open the same file for append output
        open_writers();
        // get the number of records
        size_t n_records = record_count();
        // we need to record the max value and record state during the iteration
        Key curr=0, prev=0;
        Value value;
        bool missing_records = false;
        // go through the records of the file and write records [k_i, 0x0] for each k_i that we don't see up to the max record
        for (size_t i = 0; i < n_records+1; ++i) {
            if (i == n_records) {
                // handle the max record case
                curr = max_key+1;
            } else {
                curr = read_key();
                value = read_value();
            }
            //std::cerr << "seeing " << curr << " " << value << std::endl;
            while (prev+1 < curr) {
                ++prev;
                //std::cerr << "appending " << prev << " " << nullvalue << std::endl;
                missing_records = true;
                append(prev, nullvalue);
            }
            //append(curr, value); // no need to append as we already have it!
            prev = curr;
        }
        // 
        // we have to sort again if we found any empty records
        if (missing_records) {
            sort();
        }
    }

    // index
    void index(Key new_max = 0) {
        if (new_max) max_key = new_max;
        sort();
        pad();
        open_readers();
        size_t n_records = record_count();
        sdsl::bit_vector key_bv(n_records+1);
        // record the key starts
        Key last = nullkey, curr = nullkey;
        Value val = nullvalue;
        //reader.read((char*)&last, sizeof(Key));
        //key_bv[0] = 1;
        for (size_t i = 0; i < n_records; ++i) {
            curr = read_key();
            val = read_value();
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
        auto& reader = get_reader();
        reader.seekg(n*record_size);
        return read_key();
    }

    Value nth_value(size_t n) {
        auto& reader = get_reader();
        reader.seekg(n*record_size+sizeof(Key));
        return read_value();
    }

    void for_each_pair(const std::function<void(const Key&, const Value&)>& lambda) {
        open_readers(); // open or seek to beginning
        Key key;
        Value value;
        size_t n_records = record_count();
#pragma omp parallel for
        for (size_t i = 0; i < n_records; ++i) {
            key = read_key();
            value = read_value();
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
        std::unordered_set<Value> seen;
        for_values_of(key, [&seen, &lambda](const Value& value) {
                if (!seen.count(value)) {
                    lambda(value);
                    seen.insert(value);
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
        if (!readers.size()) open_readers();
        size_t i = key_cbv_select(key);
        size_t j = key_cbv_select(key+1);
        auto& reader = get_reader();
        reader.seekg(i*record_size);
        for ( ; i < j; ++i) {
            reader.ignore(sizeof(Key));
            Value value;
            reader.read((char*)&value, sizeof(Value));
            if (!is_null(value)) {
                lambda(value);
            }
        }
    }
};

}

#endif
