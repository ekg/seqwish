# In-place Parallel Super Scalar Samplesort (IPS⁴o)

This is the implementation of the algorithm IPS⁴o presented in the paper [Engineering In-place (Shared-memory) Sorting Algorithms](https://arxiv.org/abs/2009.13569),
which contains an in-depth description of its inner workings, as well as an extensive experimental performance evaluation.
Here's the abstract:

> We present new sequential and parallel sorting algorithms that now
> represent the fastest known techniques for a wide range of input
> sizes, input distributions, data types, and machines. Somewhat
> surprisingly, part of the speed advantage is due to the additional
> feature of the algorithms to work in-place, i.e., they do not need a
> significant amount of space beyond the input array. Previously, the
> in-place feature often implied performance penalties. Our main
> algorithmic contribution is a blockwise approach to in-place data
> distribution that is provably cache-efficient.  We also parallelize
> this approach taking dynamic load balancing and memory locality into
> account.
>
> Our new comparison-based algorithm *In-place Superscalar Samplesort
> (IPS⁴o)*, combines this technique with branchless decision
> trees. By taking cases with many equal elements into account and
> by adapting the distribution degree dynamically, we obtain a
> highly robust algorithm that outperforms the best previous
> in-place parallel comparison-based sorting algorithms by almost a
> factor of three. That algorithm also outperforms the best
> comparison-based competitors regardless of whether we consider
> in-place or not in-place, parallel or sequential settings.
>
> Another surprising result is that IPS⁴o even outperforms the best
> (in-place or not in-place) integer sorting algorithms in a wide
> range of situations. In many of the remaining cases (often involving
> near-uniform input distributions, small keys, or a sequential
> setting), our new *In-place Parallel Super Scalar Radix Sort
> (IPS²Ra)* turns out to be the best algorithm.
>
> Claims to have the -- in some sense -- "best" sorting algorithm can
> be found in many papers which cannot all be true.  Therefore, we
> base our conclusions on an extensive experimental study involving a
> large part of the cross product of 21 state-of-the-art sorting
> codes, 6 data types, 10 input distributions, 4 machines, 4 memory
> allocation strategies, and input sizes varying over 7 orders of
> magnitude. This confirms the claims made about the robust
> performance of our algorithms while revealing major performance
> problems in many competitors outside the concrete set of
> measurements reported in the associated publications. This is
> particularly true for integer sorting algorithms giving one reason
> to prefer comparison-based algorithms for robust general-purpose
> sorting.

An initial version of IPS⁴o has been described in our [publication](https://drops.dagstuhl.de/opus/volltexte/2017/7854/pdf/LIPIcs-ESA-2017-9.pdf) on the 25th Annual European Symposium on Algorithms.

## Usage

Clone this repository and check out its submodule

```bash
git clone --recurse-submodules https://github.com/ips4o/ips4o.git
```

or use the following commands instead if you want to include this repository as a submodule:

```bash
git submodule add https://github.com/ips4o/ips4o.git
git submodule update --recursive --init
```

IPS⁴o provides a CMake library for simple usage:

```CMake
add_subdirectory(<path-to-the-ips4o-repository>)
target_link_libraries(<your-target> PRIVATE ips4o)
```

A minimal working example:

```C++
#include "ips4o.hpp"

// sort sequentially
ips4o::sort(begin, end[, comparator]);

// sort in parallel (uses OpenMP if available, std::thread otherwise)
ips4o::parallel::sort(begin, end[, comparator]);
```

The parallel version of IPS⁴o requires 16-byte atomic compare-and-exchange instructions to run the fastest.
Most CPUs and compilers support 16-byte compare-and-exchange instructions nowadays.
If the CPU in question does so, IPS⁴o uses 16-byte compare-and-exchange instructions when you set your CPU correctly (e.g., `-march=native`) or when you enable the instructions explicitly (`-mcx16`).
In this case, you also have to link against GCC's libatomic (`-latomic`).
Otherwise, we emulate some 16-byte compare-and-exchange instructions with locks which may slightly mitigate the performance of IPS⁴o.

If you use the CMake example shown above, we automatically optimize IPS⁴o for the native CPU (e.g., `-march=native`).
You can disable the CMake property `IPS4O_OPTIMIZE_FOR_NATIVE` to avoid native optimization and you can enable the CMake property `IPS4O_USE_MCX16` if you compile with GCC or Clang to enable 16-byte compare-and-exchange instructions explicitly.

IPS⁴o uses C++ threads if not specified otherwise.
If you prefer OpenMP threads, you need to enable OpenMP threads, e.g., enable the CMake property `IPS4O_USE_OPENMP` or add OpenMP to your target.
If you enable the CMake property `DISABLE_IPS4O_PARALLEL`, most of the parallel code will not be compiled and no parallel libraries will be linked.
Otherwise, CMake automatically enables C++ threads (e.g., `-pthread`) and links against TBB and GCC's libatomic. (Only when you compile your code for 16-byte compare-and-exchange instructions you need libatomic.)
Thus, you need the Thread Building Blocks (TBB) library to compile and execute the parallel version of IPS⁴o.
We search for TBB with `find_package(TBB REQUIRED)`.
If you want to execute IPS⁴o in parallel but your TBB library is not accessible via `find_package(TBB REQUIRED)`, you can still compile IPS⁴o with parallel support. 
Just enable the CMake property `DISABLE_IPS4O_PARALLEL`, enable C++ threads for your own target and link your own target against your TBB library (and also link your target against libatomic if you want 16-byte atomic compare-and-exchange instruction support).

If you do not set a CMake build type, we use the build type `Release` which disables debugging (e.g., `-DNDEBUG`) and enables optimizations (e.g., `-O3`).

Currently, the code does not compile on Windows.

## Licensing

IPS⁴o is free software provided under the BSD 2-Clause License described in the [LICENSE file](LICENSE). If you use this implementation of IPS⁴o in an academic setting please cite the paper [Engineering In-place (Shared-memory) Sorting Algorithms](https://arxiv.org/abs/2009.13569) using the BibTeX entry

```bibtex 
@misc{axtmann2020engineering,
  title =	 {Engineering In-place (Shared-memory) Sorting Algorithms},
  author =	 {Michael Axtmann and Sascha Witt and Daniel Ferizovic and Peter Sanders},
  howpublished = {Computing Research Repository (CoRR)},
  year =	 {Sept. 2020},
  archivePrefix ={arXiv},
  eprint =	 {2009.13569},
}
```
