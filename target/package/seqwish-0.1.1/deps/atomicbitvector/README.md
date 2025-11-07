# atomic bitvector

An atomic bitset/bitvector with size determined at runtime

## overview

This header-only library encodes a bitvector class with size fixed at runtime.
Atomic operations allow for concurrent access and modification of the bitset.
Such a structure can help with coordinating parallel processing of a given fixed set of entities, and has the advantage of only requiring O(1) bit per entry.

The `atomic_bv_t` class is a straightforward extension of [`ConcurrentBitSet`](https://github.com/facebook/folly/blob/f296a1d04b5bd13bbfcce3fa7644ec3932454d62/folly/ConcurrentBitSet.h) from Facebook's folly C++ library.
It wraps the atomic type with a class that allows them to be copied, and these wrapped atomic types are then used to build a vector whose size is determined at runtime.

## tests

Build the test using:

```
cmake -H. -Bbuild && cmake --build build -- -j4
```

And run it with:

```
n_bits=10000000
n_threads=4
time bin/test_atomicbitvector $n_bits $n_threads
```

In the test, `$n_threads` write randomly to a bitvector of `$n_bits`.
A reporting thread runs over the bitvector, checking when all bits have been set.
Parallel speedup can be confirmed by setting `$n_threads` to 1, 2, 3 ... up to the number of parallel processes allowed on the system.

## example usage

```c++
#include <atomic_bitvector.hpp>
#include <thread>

int main(int argc, char** argv) {
    atomicbitvector::atomic_bv_t x(1000);
    size_t i = 100;
    x.test(i);  // false
    x.set(i);   // false (returns previous value)
    x.test(i);  // true
    x.reset(i); // true (returns previous value)n
    x.test(i);  // false
    return 0;
}
```

All operations are atomic and threadsafe.

## license

Apache2
