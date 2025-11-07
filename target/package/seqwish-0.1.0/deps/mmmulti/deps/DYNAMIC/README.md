DYNAMIC: a succinct and compressed fully-dynamic data structures library
===============

### Contributors

The main contributors to the library are: 

- Alan Kuhnle. akuhnle418@gmail.com
- Nicola Prezza (creator).  nicola.prezza@gmail.com

The library has also received many contributions from other researchers: we wish to thank Mikhail Karasikov, Erik Garrison, and Chris Barber for adding many useful features. We are also very grateful to Adam Novak, Jan Studený, Uula Ulkuniemi, and Michael R. Crusoe for useful bug reports and suggestions. 

Please cite this library as: 

@inproceedings{prezza2017framework,  
  title={A Framework of Dynamic Data Structures for String Processing},  
  author={Prezza, Nicola},  
  booktitle={International Symposium on Experimental Algorithms},  
  year={2017},  
  organization={Leibniz International Proceedings in Informatics (LIPIcs)}  
}

### Data structures

**NEW: Delete (Remove) operations are now available on most structures! (indel = insertions + deletions)**

This library offers space- and time-efficient implementations of some basic succinct/compressed dynamic data structures. DYNAMIC features:

- A succinct Searchable Partial Sums with **Indels** (SPSI) structure representing a list of integers s_1, s_2, ..., s_m. Space: about 1.2 * m * ( log(M/m) + log log m ) bits, where M = m + s_1 + s_2 + ... + s_m. The structure supports also update operations (i.e. s_i = s_i + delta).
- A Succinct dynamic bitvector supporting rank/select/access/**Indel** (RSAI) operations. Space: about 1.2 * n bits.
- A gap-compressed dynamic bitvector supporting rank/select/access/**Indel** operations. Space: about 1.2 * b * ( log(n/b) + log log b ) bits,  b being the number of bits set and n being the bitvector length. All operations take log(b) time.
- A dynamic sparse vector (of integers) with access/**Indel** operations.
- A dynamic string supporting rank/select/access/**Indel** operations. The user can choose at construction time between fixed-length/gamma/Huffman encoding of the alphabet. All operations take log(n) * log(sigma) time (or log(n) * H0 with Huffman encoding).
- A run-length encoded dynamic string supporting rank/select/access/insert operations (removes are not yet implemented). Space: approximately R*(1.2 * log(sigma) + 2.4 * (log(n/R)+log log R) ) bits, where R is the number of runs in the string. All operations take log(R) time.
- A dynamic (left-extend only) entropy/run-length compressed BWT
- A dynamic (left-extend only) entropy/run-length compressed FM-index. This structure consists in the above BWT + a dynamic suffix array sampling

### Algorithms

- Two algorithms to build LZ77 in repetition-aware RAM working space. Both algorithms use a run-length encoded BWT with sparse Suffix array sampling. The first algorithm stores 2 SA samples per BWT run. The second algorithm (much more space efficient) stores 1 SA sample per LZ factor. From the papers "Computing LZ77 in Run-Compressed Space", Alberto Policriti and Nicola Prezza, DCC2016 and "
LZ77 Computation Based on the Run-Length Encoded BWT", Alberto Policriti and Nicola Prezza (Algorithmica)
- An algorithm to build the BWT in run-compressed space
- An algorithm to build LZ77 in nH0(2+o(1)) space and n * log n * H0 time. From the paper "Fast Online Lempel-Ziv Factorization in Compressed Space", Alberto Policriti and Nicola Prezza, SPIRE2015
- An algorithm to build the BWT in high-order compressed space. The algorithm runs in O(n * H_k * log log n) average-case time (e.g. good for DNA) and O(n * H_k * log n) worst-case time. From the paper "Average linear time and compressed space construction of the Burrows-Wheeler transform"
Policriti A., Gigante N. and Prezza N., LATA 2015 (the paper discusses a theoretically faster variant)

The SPSI structure is the building block on which all other structures are based. This structure is implemented with cache-efficient B-trees.

### TODO: 

- Implement remove operations on rle_string, bwt, and fm_index
- Dynamic wavelet matrices
- Implement a good memory allocator. At the moment the default allocator is used, which results in about 25% of memory being wasted due to fragmentation
- Geometric data structures (predecessor/2D range search)
- Batch inserts

### Download

> git clone https://github.com/nicolaprezza/dynamic

### Compile

Thre library features some example executables. To compile them, firstly create and enter a bin/ directory

> mkdir bin; cd bin

Then, launch cmake as (default build type is release):

> cmake ..

Finally, build the executables:

> make

The above command creates the executables in the bin directory. 

### Usage

The header include/dynamic/dynamic.hpp contains all type definitions and is all you need to include in your code. The folder algorithms/ contains some algorithms implemented with the library's structures. This is a snapshot of dynamic.hpp:

    /*
     * a succinct searchable partial sum with inserts implemented with cache-efficient
     * B trees.
     */
    typedef spsi<packed_vector,256,16> packed_spsi;

    /*
     * dynamic gap-encoded bitvector
     */
    typedef gap_bitvector<packed_spsi> gap_bv;

    /*
     * dynamic succinct bitvector (about 1.2n bits)
     */
    typedef succinct_bitvector<spsi<packed_vector,8192,16> > suc_bv;

    /*
     * succinct/compressed dynamic string implemented with wavelet trees.
     * user can choose (at construction time) between fixed-length / gamma / Huffman encoding of characters.
     */
    typedef wt_string<suc_bv> wt_str;

    /*
     * run-length encoded (RLE) string. This string uses 1 sparse bitvector
     * for all runs, one dynamic string for run heads, and sigma sparse bitvectors (one per character)
     */
    typedef rle_string<gap_bv, wt_str> rle_str;

    /*
     * RLE string implemented with a run-length encoded wavelet tree. Each
     * WT node is run-length encoded. 
     */
    typedef wt_string<rle_str> wtrle_str;

    /*
     * succinct/compressed BWT (see description of wt_str)
     */
    typedef bwt<wt_str,rle_str> wt_bwt;

    /*
     * run-length encoded BWT
     */
    typedef bwt<rle_str,rle_str> rle_bwt;

    /*
     * dynamic sparse vector: <= m*k + O(m log n/m) bits of space, where k is the maximum
     * number of bits of any integer > 0 and n is the total number of integers.
     */
    typedef sparse_vector<packed_spsi,gap_bv> sparse_vec;

    /*
     * dynamic succinct/entropy compressed FM index. BWT positions are
     * marked with a succinct bitvector
     *
     * roughly n*H0 + n + (n/k)*log n bits of space, where k is the SA sample  rate
     *
     */
    typedef fm_index<wt_bwt, suc_bv, packed_spsi> wt_fmi;

    /*
     * dynamic run-length encoded FM index. BWT positions are
     * marked with a gap-encoded bitvector.
     *
     * roughly 2.4*R*log(n/R) + 2.4 log log R + 1.2*R*log(sigma) + (n/k)*log n bits of  space, where
     * k is the SA sample rate and R is the number of runs in the BWT
     *
     */
    typedef fm_index<rle_bwt, gap_bv, packed_spsi> rle_fmi;
