# mmmulti

memory-mapped multimap, multiset, and (implicit) interval tree

## rationale

Sometimes you have a lot of plain-old data, but you need random access to it.
These header-only classes combine memory-mapped files with high-performance parallel sorting and appropriate indexing strategies to support very large (>memory but <disk) multimaps, multisets, and interval trees.

## mmmulti::map memory-mapped multimap

This implements a memory backed multimap intended for use where:

- your keys are integers, or can be mapped to dense range of integers,
- the memory mapped file is on fast storage, like an SSD (although this is not a requirement),
- you have arbitrary values of fixed size (e.g. structs, other POD types) that can be sorted,
- you don't need dynamic updates of the table,
- and you are likely to run out of memory of you use a traditional map or hash table,
- but you can handle approximately 1 bit per record in RAM.

These may seem to be very specific, but many problems can be mapped into a dense integer set.
`mmmulti::map` developed first as a data structure to support [seqwish](https://github.com/ekg/seqwish), which uses it to generate precise variation graphs from pairwise alignments between collections of sequences.
As this multimap forms a key data processing kernel in the algorithm, it can scale to extremely large problem sizes, limited only by available disk space.
Although performance is much slower than an in-memory structure, we are virtually guaranteed to be able to complete the compute.

### usage

To construct the `mmmulti::map`:

```c++
#include "mmmultimap.hpp"

mmmulti::map<uint64_t, uint64_t> mmap("temp.dat");
```

We can then add key pairs:

```c++
mmap.open_writer(); // required before adding keys, opens a writer/coordinator thread
mmap.append(key1, value1);
mmap.append(key1, value2);
mmap.append(key1, value3);
```

Calls to `append` are threadsafe once `open_writer` has been called.

To query the `mmmulti::map`, first index it, providing the maximum key to expect (remember, we're working on dense keys!):

```c++
mmap.index(num_threads, max_key);
```

If we index without providing a maximum key, like this:

```c++
mmap.index(num_threads);
```

... then we don't pad the multimap, and we can only enumerate the keys inside using e.g. `for_each_pair`.
This can have some advantages, such as not requiring that our values be non-null (positive integers for instance, or structs with non-null entries).
This allows us to simply use the sorted array functionality of the `mmmulti::map`.

#### indexing algorithm

If `max_key` is specified, we first pads the records with one key/null record per key in the range of `[0,max_key]`.
Then, we memory map this file and apply the [ips4o](https://github.com/SaschaWitt/ips4o) in-place parallel super scalar samplesort to the memory mapped buffer to order the records.
When padded (`max_key` is specified), the index is completed by building a bitvector of length `max_key`, marking 1 for the first instance of a key in the sorted order, and building an auxiliary data structure to support `select_1(n)` queries on it that allow us to find the records associated with a given key.
Without this index, we can only iterate through the key/value pairs.

#### supported queries

It is now possible to iterate through the records in order with `for_each_pair`.
We can look up the nth key or value with `nth_key` and `nth_value`, which can be used to enable parallel traversal of the keys externally.
And we can iterate over the values or unique values of a given key with `for_values_of` and `for_unique_values_of`.

## mmmulti::set memory-mapped multiset with iteration

This is similar to `mmmulti::map`, but useful where random access to values is not required, and where contiguity of keys is not possible.
Unlike the `mmmulti::map`, arbitrary fixed-length data types are allowed, but they must have `operator<` defined to enable sorting.
It drops the index structures and padding, saving space, but preserves the same API.
The `mmmulti::set` only provides iteration across its key space.
As such, it's useful where we need to collect and count a set of entities.
Random access is not currently supported (TODO: it would be easy to implement using binary search, but current applications do not require it).

### usage

To construct the `mmmulti::set`:

```c++
#include "mmmultiset.hpp"

mmmulti::set<uint64_t> mset("temp.dat");
```

We can then add values:

```c++
mset.open_writer(); // required before adding keys, opens a writer/coordinator thread
mset.append(value1);
mset.append(value2);
mset.append(value3);
```

Calls to `append` are threadsafe once `open_writer` has been called.

To use the `mmmulti::set`, first index it:

```c++
mset.index(num_threads);
```

#### indexing algorithm

Indexing closes the writer, memory maps the backing file and applies the [ips4o](https://github.com/SaschaWitt/ips4o) in-place parallel super scalar samplesort to the memory mapped buffer to order the records.

#### supported queries

It is now possible to iterate through the records in order with `for_each_value`, along with their counts with `for_each_value_count`, or just unique values with `for_each_unique_value`.

## mmmulti::iitree memory-mapped implicit interval tree

The implicint interval tree data structure sorts a collection of intervals into a linear array ([cgranges](https://github.com/lh3/cgranges)).
Tree traversal is achieved by jumping between array indexes.
`mmmulti::iitree` implements this model on top of a disk-backed memory-mapped array, and uses the [ips4o](https://github.com/SaschaWitt/ips4o) in-place parallel super scalar to speed up sorting and index generation.
Usage is similar to other classes in `mmmulti`.

### usage

To construct the `mmmulti::set`:

```c++
#include "mmiitree.hpp"

# template arguments are range integer type and stored value
mmmulti::iitree<uint64_t, Data> tree("temp.dat");
```

We can then add ranges:

```c++
tree.open_writer(); // required before adding intervals, opens a writer/coordinator thread
tree.add(start1, end1, data1);
tree.add(start2, end2, data2);
tree.add(start3, end3, data3);
```

Calls to `add` are threadsafe once `open_writer` has been called.

To use the `mmmulti::iitree`, first index it:

```c++
tree.index(num_threads);
```

#### indexing algorithm

Indexing closes the writer, memory maps the backing file and applies the [ips4o](https://github.com/SaschaWitt/ips4o) in-place parallel super scalar samplesort to the memory mapped buffer to order the records.
The indexing procedure from [cgranges](https://github.com/lh3/cgranges) is then applied to set up implicit interval tree.

#### supported queries

To find overlaps for a given query, use `mmmulti::iitree::overlap`.
This returns a vector of range ranks in the sorted set of ranges.
We can then look up the start, end, and data fields of these in the backing tree.
For efficiency, this is done by callback.

```c++
tree.overlap(
    n, m,
    [&](const uint64_t& start,
        const uint64_t& end,
        const Data& data) {
        // process record
        std::cout << "start " << start << std::endl;
        std::cout << "end   " << end   << std::endl;
        std::cout << "data  " << data  << std::endl;
    });
```

It's also possible to iterate through the ranges with `for_each_entry` and also with their counts, where duplicates are present with `for_each_entry_count`.

## building and testing

`mmmulti`'s classes are intended to be used as libraries in other applications.
The namespace can be included easily as a CMake ExternalProject.
A test utility is included to demonstrate usage.
To build it:

```
cmake -H. -Bbuild && cmake --build build -- -j 4
```

And to run some tests:

```
bin/mmmulti -t x -s 10000000 -t 4
10040123 keys
15099011 values
10688272 unique pairs
rm -f x # removes test file
```

## development

By adding a PMHF to the frontend, it should be possible to project arbitrary key sets into the dense range required by `mmmulti::map` with only a few bits of overhead per entry.
This would also obviate the need to pad our key space.

## author

Erik Garrison <erik.garrison@gmail.com>

## license

MIT
