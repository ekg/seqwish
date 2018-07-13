#ifndef SEQINDEX_HPP_INCLUDED
#define SEQINDEX_HPP_INCLUDED

#include <iostream>
#include <cstdio>
#include <string>
#include "sdsl/bit_vectors.hpp"
#include "sdsl/csa_wt.hpp"
#include "sdsl/suffix_arrays.hpp"
#include "sdsl/dac_vector.hpp"
#include "gzstream.h"

namespace seqwish {

class SeqIndex {

private:

    std::string basefilename;
    std::string seqfile;
    std::string seqnamefile;
    std::string seqidxfile;
    sdsl::dac_vector<> seq_offset_civ;   // sequence offsets (for offset and length)
    sdsl::csa_wt<> seq_name_csa;         // seq name compressed suffix array
    sdsl::sd_vector<> seq_name_cbv;      // path name starts
    sdsl::sd_vector<>::rank_1_type seq_name_cbv_rank;
    sdsl::sd_vector<>::select_1_type seq_name_cbv_select;
    uint32_t OUTPUT_VERSION = 1; // update as we change our format

public:

    SeqIndex(void) { }

    ~SeqIndex(void) { }

    // load a FASTA or FASTQ file into a file with a name index mapping name -> offset and indexed with a CSA
    // provide queries over this index that let us extract particular positions and subsequences
    void set_base_filename(const std::string& filename) {
        basefilename = filename;
        seqfile = basefilename + ".seq";
        seqnamefile = basefilename + ".seqnames";
        seqidxfile = basefilename + ".seqidx";
    }

    void build_index(const std::string& filename) {
        set_base_filename(filename);
        // read the file
        igzstream in(filename.c_str());
        std::ofstream seqnames(seqnamefile.c_str());
        std::ofstream seqout(seqfile.c_str());
        std::vector<uint64_t> seqname_offset;
        std::vector<uint64_t> seq_offset;
        bool input_is_fasta=false, input_is_fastq=false;
        // look at the first character to determine if it's fastq or fasta
        std::string line;
        std::getline(in, line);
        if (line[0] == '>') {
            input_is_fasta = true;
        } else if (line[0] == '@') {
            input_is_fastq = true;
        } else {
            std::cerr << "unknown file format given to SeqIndex" << std::endl;
            assert(false);
        }
        size_t seq_bytes_written = 0;
        size_t seq_names_bytes_written = 0;
        while (in.good()) {
            seqname_offset.push_back(seq_names_bytes_written);
            seq_offset.push_back(seq_bytes_written);
            line[0] = '>';
            seqnames << line;
            seq_names_bytes_written += line.size();
            std::string seq;
            // get the sequence
            if (input_is_fasta) {
                while (std::getline(in, line)) {
                    if (line[0] == '>') {
                        // seek back, this is the next sequence
                        //in.seekg(in.tellg()-(std::streamoff)(line.size()+1));
                        break;
                    } else {
                        seq.append(line);
                    }
                }
            } else if (input_is_fastq) {
                std::getline(in, seq); // sequence
                std::getline(in, line); // delimiter
                std::getline(in, line); // quality
                std::getline(in, line);
            }
            seqout << seq;
            // record where the sequence starts
            seq_bytes_written += seq.size();
        }
        in.close();
        // add the last value so we can get sequence length for the last sequence and name
        seq_offset.push_back(seq_bytes_written);
        seqname_offset.push_back(seq_names_bytes_written);
        seqnames.close();
        seqout.close();
        // mark the seq name starts vector
        sdsl::bit_vector seq_name_starts(seqname_offset.back());
        for (size_t i = 0; i < seqname_offset.size()-1; ++i) {
            seq_name_starts[i] = 1;
        }
        // build the name index
        construct(seq_name_csa, seqnamefile, 1);
        // destroy the file
        std::remove(seqnamefile.c_str());
        // build the rest of the index
        sdsl::util::assign(seq_name_cbv, sdsl::sd_vector<>(seq_name_starts));
        sdsl::util::assign(seq_name_cbv_rank, sdsl::sd_vector<>::rank_1_type(&seq_name_cbv));
        sdsl::util::assign(seq_name_cbv_select, sdsl::sd_vector<>::select_1_type(&seq_name_cbv));
        sdsl::util::assign(seq_offset_civ, sdsl::dac_vector<>(seq_offset));
    }

    size_t save(sdsl::structure_tree_node* s = NULL, std::string name = "") {
        assert(seq_name_csa.size() && seq_name_cbv.size() && seq_offset_civ.size());
        sdsl::structure_tree_node* child = sdsl::structure_tree::add_child(s, name, sdsl::util::class_name(*this));
        // open the sdsl index
        std::ofstream out(seqidxfile.c_str());
        size_t written = 0;
        out << "seqidx"; written += 9;
        uint32_t version_buffer = OUTPUT_VERSION;
        out.write((char*) &version_buffer, sizeof(version_buffer));
        written += seq_name_csa.serialize(out, child, "seq_name_csa");
        written += seq_name_cbv.serialize(out, child, "seq_name_cbv");
        written += seq_name_cbv_rank.serialize(out, child, "seq_name_cbv_rank");
        written += seq_name_cbv_select.serialize(out, child, "seq_name_cbv_select");
        written += seq_offset_civ.serialize(out, child, "seq_offset_civ");
        out.close();
        return written;
    }

    void load(const std::string& filename) {
        set_base_filename(filename);
        std::ifstream in(seqidxfile.c_str());
        std::string magic;
        in.read((char*)magic.c_str(), 6);
        uint32_t version;
        in.read((char*) &version, sizeof(version));
        assert(version == OUTPUT_VERSION);
        seq_name_csa.load(in);
        seq_name_cbv.load(in);
        seq_name_cbv_rank.load(in);
        seq_name_cbv_select.load(in);
        seq_offset_civ.load(in);
        in.close();
    }
};

}

#endif
