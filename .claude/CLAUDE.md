# Seqwish Project Requirements

## CRITICAL: Big O Complexity Preservation

**NON-NEGOTIABLE REQUIREMENT:**

When migrating C++ code to Rust, you **MUST** use algorithms and data structures with the **EXACT SAME Big O bounds** in both time and space as the original C++ implementation.

### Rules

1. **No space complexity trade-offs**
   - If C++ uses a compressed suffix array (CSA), Rust MUST use CSA or equivalent
   - If C++ uses succinct bitvectors, Rust MUST use succinct bitvectors
   - If C++ uses O(n) space, Rust MUST use O(n) space
   - "Practical" arguments (e.g., "only 100KB extra for small datasets") are NOT acceptable

2. **No time complexity trade-offs**
   - If C++ has O(log n) lookup, Rust MUST have O(log n) lookup
   - If C++ has O(n log n) sorting, Rust MUST have O(n log n) sorting
   - Constant factors matter for performance, but Big O bounds are non-negotiable

3. **Research first**
   - Before implementing, find Rust crates that match C++ Big O bounds
   - If no Rust equivalent exists, use FFI to wrap C++ implementation
   - Document why a particular Rust crate matches the complexity requirements

4. **Never compromise**
   - There is NO case where it's acceptable to change Big O bounds
   - If stuck, ask the user before proceeding
   - If uncertain about complexity, research or ask

### Examples of UNACCEPTABLE Changes

❌ **WRONG:** Replace CSA (O(m log n) space) with HashMap (O(n) space)
✅ **CORRECT:** Use Rust CSA implementation or FFI to sdsl::csa_wt

❌ **WRONG:** Replace succinct bitvector (O(n) bits) with Vec<u64> (O(n) words)
✅ **CORRECT:** Use vers-vecs or simple-sds with equivalent rank/select support

❌ **WRONG:** Replace implicit interval tree (O(n) space) with pointer-based tree (O(n log n) space)
✅ **CORRECT:** Port the implicit tree structure or use equivalent implementation

### Why This Matters

Seqwish is designed for pangenomes with:
- Gigabytes to terabytes of sequence data
- Millions to billions of alignment intervals
- Memory-constrained cluster environments

Even "small" constant factor increases in space usage can make the difference between:
- Running on available hardware vs requiring bigger machines
- Processing a dataset vs running out of memory
- Being usable vs being abandoned

**This is not negotiable.**
