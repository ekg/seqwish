# paryfor

a parallel_for implementation based on atomic queues

## usage

The `parallel_for` templates replace `#pragma omp parallel for`.
We use a callback that takes the iteration number.
To avoid compiler confusion, we can specify the iteration type in the template.

```c++
#include "paryfor.hpp"

int main(int argc, char** argv) {
    uint64_t todo_count = std::stoul(argv[1]);
    int thread_count = std::stoi(argv[2]);
    int chunk_size = std::stoi(argv[3]);
    paryfor::parallel_for<uint64_t>(
        0, todo_count, thread_count, chunk_size,
        [&](uint64_t i, int tid) {
            // do work
        });
}
```

We don't need to use `chunk_size`.
If omitted, we use a template that queues single iteration values, rather than ranges.

We can also pass a function that does not take its thread id as an argument.

## thanks

This utility depends on Maxim Egorushkin's atomic_queue library, which is included in the single header file `paryfor.hpp`.

## author

Erik Garrison

## license

MIT
