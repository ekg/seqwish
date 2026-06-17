//! # seqwish - A variation graph inducer
//!
//! Seqwish builds variation graphs from pairwise sequence alignments.
//! It transforms a collection of sequences and their all-to-all alignments
//! into a graph representation that captures the variation between the sequences.
//!
//! ## Overview
//!
//! The algorithm proceeds in several stages:
//!
//! 1. **Sequence Indexing** - Load and index input sequences
//! 2. **Alignment Processing** - Parse and index PAF alignments
//! 3. **Transitive Closure** - Compute equivalence classes of aligned positions
//! 4. **Node Compaction** - Merge non-bifurcating regions into single nodes
//! 5. **Link Derivation** - Extract edges between nodes
//! 6. **GFA Emission** - Output the variation graph in GFA format
//!
//! ## Example
//!
//! ```rust,no_run
//! use seqwish::seqindex::SeqIndex;
//! use std::sync::{Arc, Mutex};
//!
//! // Build a sequence index
//! let mut seqidx = SeqIndex::new();
//! seqidx.build_index("sequences.fa").unwrap();
//! ```
//!
//! ## Command-line Usage
//!
//! ```bash
//! seqwish -s sequences.fa -p alignments.paf -g output.gfa
//! ```
//!
//! ## Features
//!
//! - Memory-safe parallel processing
//! - Disk-backed data structures for scalability
//! - Produces GFA v1.0 format output
//! - Compatible with standard pangenome tools

// Allow clippy warnings for FFI functions that dereference raw pointers
// These are inherently unsafe operations but are needed for C compatibility
#![allow(clippy::not_unsafe_ptr_arg_deref)]
// Allow dead code - some functions are exported for FFI or future use
#![allow(dead_code)]
// Allow complex types and many arguments in certain contexts
#![allow(clippy::too_many_arguments)]
#![allow(clippy::type_complexity)]

use std::ffi::{c_char, CStr, CString};
use std::ptr;
use std::sync::Arc;

use bitvec::prelude::*;

pub mod alignments;
pub mod cigar;
pub mod compact;
pub mod dna;
pub mod dset64;
pub mod dset64_asm;
pub mod dset64_unsafe;
pub mod gfa;

#[cfg(feature = "cuda")]
pub mod gpu;

#[cfg(not(feature = "cuda"))]
pub mod gpu {
    pub fn is_cuda_available() -> bool {
        false
    }
    pub fn gpu_union_find(
        _num_elements: usize,
        _edges: &[(u32, u32)],
        _verbose: bool,
    ) -> Option<Vec<u32>> {
        None
    }
    pub struct GpuRunner;
    impl GpuRunner {
        pub fn new() -> Option<Self> {
            None
        }
        pub fn gpu_union_find(
            &mut self,
            _n: usize,
            _e: &[(u32, u32)],
            _v: bool,
        ) -> Option<Vec<u32>> {
            None
        }
    }
}

pub mod intervaltree;
pub mod links;
pub mod mmap;
pub mod paf;
pub mod pos;
pub mod seqindex;
pub mod sxs;
pub mod tempfile;
pub mod time;
pub mod transclosure;
pub mod utils;
pub mod version;

/// Returns the version string of the Rust component
#[no_mangle]
pub extern "C" fn seqwish_rust_version() -> *const c_char {
    b"0.1.0-rust\0".as_ptr() as *const c_char
}

/// Simple test function to verify FFI is working
#[no_mangle]
pub extern "C" fn seqwish_rust_add(a: i32, b: i32) -> i32 {
    a + b
}

// FFI wrappers for tempfile module

/// Create a temporary file. Returns a C string that must be freed with temp_file_free_string.
/// Returns NULL on error.
#[no_mangle]
pub extern "C" fn temp_file_create(base: *const c_char, suffix: *const c_char) -> *mut c_char {
    if base.is_null() || suffix.is_null() {
        return ptr::null_mut();
    }

    let base_str = unsafe {
        match CStr::from_ptr(base).to_str() {
            Ok(s) => s,
            Err(_) => return ptr::null_mut(),
        }
    };

    let suffix_str = unsafe {
        match CStr::from_ptr(suffix).to_str() {
            Ok(s) => s,
            Err(_) => return ptr::null_mut(),
        }
    };

    match tempfile::create(base_str, suffix_str) {
        Ok(path) => {
            let path_str = path.to_string_lossy().to_string();
            match CString::new(path_str) {
                Ok(c_string) => c_string.into_raw(),
                Err(_) => ptr::null_mut(),
            }
        }
        Err(_) => ptr::null_mut(),
    }
}

/// Remove a temporary file
#[no_mangle]
pub extern "C" fn temp_file_remove(filename: *const c_char) {
    if filename.is_null() {
        return;
    }

    let filename_str = unsafe {
        match CStr::from_ptr(filename).to_str() {
            Ok(s) => s,
            Err(_) => return,
        }
    };

    tempfile::remove(std::path::Path::new(filename_str));
}

/// Set temp directory
#[no_mangle]
pub extern "C" fn temp_file_set_dir(dir: *const c_char) {
    if dir.is_null() {
        return;
    }

    let dir_str = unsafe {
        match CStr::from_ptr(dir).to_str() {
            Ok(s) => s,
            Err(_) => return,
        }
    };

    tempfile::set_dir(dir_str);
}

/// Get temp directory. Returns a C string that must be freed with temp_file_free_string.
#[no_mangle]
pub extern "C" fn temp_file_get_dir() -> *mut c_char {
    let dir = tempfile::get_dir();
    let dir_str = dir.to_string_lossy().to_string();

    match CString::new(dir_str) {
        Ok(c_string) => c_string.into_raw(),
        Err(_) => ptr::null_mut(),
    }
}

/// Set whether to keep temp files
#[no_mangle]
pub extern "C" fn temp_file_set_keep_temp(setting: bool) {
    tempfile::set_keep_temp(setting);
}

/// Clean up all temp files and parent directory
/// Call this after each graph build when processing multiple graphs in the same process
#[no_mangle]
pub extern "C" fn temp_file_cleanup() {
    tempfile::cleanup();
}

/// Free a string returned by temp_file functions
#[no_mangle]
pub extern "C" fn temp_file_free_string(s: *mut c_char) {
    if !s.is_null() {
        unsafe {
            let _ = CString::from_raw(s);
        }
    }
}

// FFI wrappers for pos module

/// Create a position from offset and orientation
#[no_mangle]
pub extern "C" fn pos_make_pos_t(offset: u64, is_rev: bool) -> u64 {
    pos::make_pos_t(offset, is_rev)
}

/// Extract offset from position
#[no_mangle]
pub extern "C" fn pos_offset(pos: u64) -> u64 {
    pos::offset(pos)
}

/// Check if position is reverse
#[no_mangle]
pub extern "C" fn pos_is_rev(pos: u64) -> bool {
    pos::is_rev(pos)
}

/// Increment position
#[no_mangle]
pub extern "C" fn pos_incr_pos(pos: *mut u64) {
    if !pos.is_null() {
        unsafe {
            pos::incr_pos(&mut *pos);
        }
    }
}

/// Increment position by N
#[no_mangle]
pub extern "C" fn pos_incr_pos_by(pos: *mut u64, by: usize) {
    if !pos.is_null() {
        unsafe {
            pos::incr_pos_by(&mut *pos, by);
        }
    }
}

/// Decrement position
#[no_mangle]
pub extern "C" fn pos_decr_pos(pos: *mut u64) {
    if !pos.is_null() {
        unsafe {
            pos::decr_pos(&mut *pos);
        }
    }
}

/// Decrement position by N
#[no_mangle]
pub extern "C" fn pos_decr_pos_by(pos: *mut u64, by: usize) {
    if !pos.is_null() {
        unsafe {
            pos::decr_pos_by(&mut *pos, by);
        }
    }
}

/// Reverse position orientation
#[no_mangle]
pub extern "C" fn pos_rev_pos_t(pos: u64) -> u64 {
    pos::rev_pos_t(pos)
}

/// Convert position to string (returns C string that must be freed)
#[no_mangle]
pub extern "C" fn pos_to_string_c(pos: u64) -> *mut c_char {
    let s = pos::pos_to_string(pos);
    match CString::new(s) {
        Ok(c_string) => c_string.into_raw(),
        Err(_) => ptr::null_mut(),
    }
}

// FFI wrappers for dna module

/// Get complement of a single DNA base
#[no_mangle]
pub extern "C" fn dna_complement(c: u8) -> u8 {
    dna::complement(c)
}

/// Reverse complement a DNA sequence (allocates new string that must be freed)
#[no_mangle]
pub extern "C" fn dna_reverse_complement(seq: *const c_char, len: usize, out: *mut c_char) {
    if seq.is_null() || out.is_null() {
        return;
    }

    unsafe {
        let slice = std::slice::from_raw_parts(seq as *const u8, len);
        let rc = dna::reverse_complement(slice);
        std::ptr::copy_nonoverlapping(rc.as_ptr(), out as *mut u8, len);
    }
}

/// Reverse complement a DNA sequence in place
#[no_mangle]
pub extern "C" fn dna_reverse_complement_in_place(seq: *mut c_char, len: usize) {
    if seq.is_null() {
        return;
    }

    unsafe {
        let slice = std::slice::from_raw_parts_mut(seq as *mut u8, len);
        dna::reverse_complement_in_place(slice);
    }
}

// FFI wrappers for cigar module

/// Opaque handle to CIGAR vector
pub struct CigarHandle {
    cigar: Vec<cigar::CigarOp>,
}

/// Opaque handle to SeqIndex
pub struct SeqIndexHandle {
    seqidx: Arc<seqindex::SeqIndex>,
}

/// Opaque handle to IITree (for node/path iitrees that use RwLock)
pub struct IITreeHandle {
    iitree: Arc<std::sync::RwLock<crate::intervaltree::AdaptiveTree<u64, pos::PosT>>>,
}

/// Opaque handle to Alignment IITree (uses Mutex for writing)
pub struct AlnIITreeHandle {
    iitree: Arc<std::sync::Mutex<crate::intervaltree::AdaptiveTree<u64, pos::PosT>>>,
}

/// Parse CIGAR string and return handle to CIGAR vector
/// Returns NULL on error. Must be freed with seqwish_cigar_free.
#[no_mangle]
pub extern "C" fn cigar_from_string(s: *const c_char) -> *mut CigarHandle {
    if s.is_null() {
        return ptr::null_mut();
    }

    let s_str = unsafe {
        match CStr::from_ptr(s).to_str() {
            Ok(s) => s,
            Err(_) => return ptr::null_mut(),
        }
    };

    let cigar = cigar::cigar_from_string(s_str);
    Box::into_raw(Box::new(CigarHandle { cigar }))
}

/// Convert CIGAR vector to string
/// Returns C string that must be freed with temp_file_free_string
#[no_mangle]
pub extern "C" fn cigar_to_string(handle: *const CigarHandle) -> *mut c_char {
    if handle.is_null() {
        return ptr::null_mut();
    }

    let cigar_handle = unsafe { &*handle };
    let s = cigar::cigar_to_string(&cigar_handle.cigar);

    match CString::new(s) {
        Ok(c_string) => c_string.into_raw(),
        Err(_) => ptr::null_mut(),
    }
}

/// Get number of operations in CIGAR
#[no_mangle]
pub extern "C" fn cigar_length(handle: *const CigarHandle) -> usize {
    if handle.is_null() {
        return 0;
    }
    let cigar_handle = unsafe { &*handle };
    cigar_handle.cigar.len()
}

/// Get operation at index
/// Returns false if index out of bounds
#[no_mangle]
pub extern "C" fn cigar_get_op(
    handle: *const CigarHandle,
    index: usize,
    len_out: *mut u64,
    op_out: *mut u8,
) -> bool {
    if handle.is_null() || len_out.is_null() || op_out.is_null() {
        return false;
    }

    let cigar_handle = unsafe { &*handle };
    if index >= cigar_handle.cigar.len() {
        return false;
    }

    unsafe {
        *len_out = cigar_handle.cigar[index].len;
        *op_out = cigar_handle.cigar[index].op;
    }
    true
}

/// Free CIGAR handle
#[no_mangle]
pub extern "C" fn seqwish_cigar_free(handle: *mut CigarHandle) {
    if !handle.is_null() {
        unsafe {
            let _ = Box::from_raw(handle);
        }
    }
}

// FFI wrappers for mmap module

/// Open a file and memory-map it
/// Returns the file size on success, 0 on error
/// The buffer pointer and file descriptor are written to the provided pointers
#[no_mangle]
pub extern "C" fn mmap_open_rust(
    filename: *const c_char,
    buf_out: *mut *mut c_char,
    fd_out: *mut i32,
) -> usize {
    if filename.is_null() || buf_out.is_null() || fd_out.is_null() {
        return 0;
    }

    let filename_str = unsafe {
        match CStr::from_ptr(filename).to_str() {
            Ok(s) => s,
            Err(_) => return 0,
        }
    };

    match mmap::mmap_open(filename_str) {
        Ok(handle) => {
            unsafe {
                *buf_out = handle.ptr;
                *fd_out = handle.fd;
            }
            let size = handle.size;
            // Prevent Drop from running - we're transferring ownership to C++
            std::mem::forget(handle);
            size
        }
        Err(_) => 0,
    }
}

/// Close a memory-mapped file
#[no_mangle]
pub extern "C" fn mmap_close_rust(buf: *mut c_char, fd: i32, size: usize) {
    if buf.is_null() {
        return;
    }

    let mut handle = mmap::MmapHandle { ptr: buf, fd, size };

    mmap::mmap_close(&mut handle);
}

// FFI wrappers for utils module

/// Check if a file exists
#[no_mangle]
pub extern "C" fn file_exists(filename: *const c_char) -> bool {
    if filename.is_null() {
        return false;
    }

    let filename_str = unsafe {
        match CStr::from_ptr(filename).to_str() {
            Ok(s) => s,
            Err(_) => return false,
        }
    };

    utils::file_exists(filename_str)
}

/// Parse a number with optional suffix (k, m, g)
#[no_mangle]
pub extern "C" fn handy_parameter(value: *const c_char, default_value: f64) -> f64 {
    if value.is_null() {
        return default_value;
    }

    let value_str = unsafe {
        match CStr::from_ptr(value).to_str() {
            Ok(s) => s,
            Err(_) => return default_value,
        }
    };

    utils::handy_parameter(value_str, default_value)
}

// FFI wrappers for time module

/// Get milliseconds since Unix epoch
#[no_mangle]
pub extern "C" fn time_since_epoch_ms() -> u64 {
    time::time_since_epoch_ms()
}

// FFI wrappers for paf module

/// Parse PAF spec string, calling callback for each (filename, weight) pair
/// Callback signature: void callback(void* user_data, const char* filename, uint64_t weight)
#[no_mangle]
pub extern "C" fn parse_paf_spec(
    spec: *const c_char,
    user_data: *mut std::ffi::c_void,
    callback: Option<extern "C" fn(*mut std::ffi::c_void, *const c_char, u64)>,
) {
    if spec.is_null() || callback.is_none() {
        return;
    }

    let spec_str = unsafe {
        match CStr::from_ptr(spec).to_str() {
            Ok(s) => s,
            Err(_) => return,
        }
    };

    let callback_fn = callback.unwrap();
    for (filename, weight) in paf::parse_paf_spec(spec_str) {
        if let Ok(c_filename) = CString::new(filename) {
            callback_fn(user_data, c_filename.as_ptr(), weight);
        }
    }
}

// FFI wrappers for alignments module

/// Hash function for match parameters
#[no_mangle]
pub extern "C" fn match_hash(q: u64, t: u64, l: u64) -> u64 {
    alignments::match_hash(q, t, l)
}

/// Determine if a match should be kept based on sparsification factor
#[no_mangle]
pub extern "C" fn keep_sparse(q: u64, t: u64, l: u64, f: f32) -> bool {
    alignments::keep_sparse(q, t, l, f)
}

// FFI wrappers for paf module

/// Opaque handle to a parsed PAF row
pub struct PafRowHandle {
    row: paf::PafRow,
}

// FFI wrappers for sxs module

/// Opaque handle to a parsed SXS alignment
pub struct SxsHandle {
    aln: sxs::SxsAlignment,
}

/// Parse a PAF row from a C string line
/// Returns NULL if parsing fails
#[no_mangle]
pub extern "C" fn paf_row_parse(line: *const c_char) -> *mut PafRowHandle {
    if line.is_null() {
        return ptr::null_mut();
    }

    let line_str = unsafe {
        match CStr::from_ptr(line).to_str() {
            Ok(s) => s,
            Err(_) => return ptr::null_mut(),
        }
    };

    match paf::PafRow::from_line(line_str) {
        Some(row) => Box::into_raw(Box::new(PafRowHandle { row })),
        None => ptr::null_mut(),
    }
}

/// Free a PAF row handle
#[no_mangle]
pub extern "C" fn paf_row_free(handle: *mut PafRowHandle) {
    if !handle.is_null() {
        unsafe {
            let _ = Box::from_raw(handle);
        }
    }
}

// Field accessors for PAF row
#[no_mangle]
pub extern "C" fn paf_row_query_sequence_name(handle: *const PafRowHandle) -> *mut c_char {
    if handle.is_null() {
        return ptr::null_mut();
    }
    let row = unsafe { &(*handle).row };
    match CString::new(row.query_sequence_name.clone()) {
        Ok(s) => s.into_raw(),
        Err(_) => ptr::null_mut(),
    }
}

#[no_mangle]
pub extern "C" fn paf_row_target_sequence_name(handle: *const PafRowHandle) -> *mut c_char {
    if handle.is_null() {
        return ptr::null_mut();
    }
    let row = unsafe { &(*handle).row };
    match CString::new(row.target_sequence_name.clone()) {
        Ok(s) => s.into_raw(),
        Err(_) => ptr::null_mut(),
    }
}

#[no_mangle]
pub extern "C" fn paf_row_query_sequence_length(handle: *const PafRowHandle) -> u64 {
    if handle.is_null() {
        return 0;
    }
    unsafe { (*handle).row.query_sequence_length }
}

#[no_mangle]
pub extern "C" fn paf_row_query_start(handle: *const PafRowHandle) -> u64 {
    if handle.is_null() {
        return 0;
    }
    unsafe { (*handle).row.query_start }
}

#[no_mangle]
pub extern "C" fn paf_row_query_end(handle: *const PafRowHandle) -> u64 {
    if handle.is_null() {
        return 0;
    }
    unsafe { (*handle).row.query_end }
}

#[no_mangle]
pub extern "C" fn paf_row_query_target_same_strand(handle: *const PafRowHandle) -> bool {
    if handle.is_null() {
        return false;
    }
    unsafe { (*handle).row.query_target_same_strand }
}

#[no_mangle]
pub extern "C" fn paf_row_target_sequence_length(handle: *const PafRowHandle) -> u64 {
    if handle.is_null() {
        return 0;
    }
    unsafe { (*handle).row.target_sequence_length }
}

#[no_mangle]
pub extern "C" fn paf_row_target_start(handle: *const PafRowHandle) -> u64 {
    if handle.is_null() {
        return 0;
    }
    unsafe { (*handle).row.target_start }
}

#[no_mangle]
pub extern "C" fn paf_row_target_end(handle: *const PafRowHandle) -> u64 {
    if handle.is_null() {
        return 0;
    }
    unsafe { (*handle).row.target_end }
}

#[no_mangle]
pub extern "C" fn paf_row_num_matches(handle: *const PafRowHandle) -> u64 {
    if handle.is_null() {
        return 0;
    }
    unsafe { (*handle).row.num_matches }
}

#[no_mangle]
pub extern "C" fn paf_row_alignment_block_length(handle: *const PafRowHandle) -> u64 {
    if handle.is_null() {
        return 0;
    }
    unsafe { (*handle).row.alignment_block_length }
}

#[no_mangle]
pub extern "C" fn paf_row_mapping_quality(handle: *const PafRowHandle) -> u16 {
    if handle.is_null() {
        return 0;
    }
    unsafe { (*handle).row.mapping_quality }
}

#[no_mangle]
pub extern "C" fn paf_row_cigar(handle: *const PafRowHandle) -> *mut CigarHandle {
    if handle.is_null() {
        return ptr::null_mut();
    }
    let row = unsafe { &(*handle).row };
    Box::into_raw(Box::new(CigarHandle {
        cigar: row.cigar.clone(),
    }))
}

/// Create a new empty SXS alignment
#[no_mangle]
pub extern "C" fn sxs_new() -> *mut SxsHandle {
    Box::into_raw(Box::new(SxsHandle {
        aln: sxs::SxsAlignment::new(),
    }))
}

/// Parse SXS alignment from array of C strings (lines)
/// Returns NULL if parsing fails
#[no_mangle]
pub extern "C" fn sxs_parse_lines(lines: *const *const c_char, num_lines: usize) -> *mut SxsHandle {
    if lines.is_null() {
        return ptr::null_mut();
    }

    let mut line_vec = Vec::new();
    for i in 0..num_lines {
        unsafe {
            let line_ptr = *lines.add(i);
            if line_ptr.is_null() {
                continue;
            }
            match CStr::from_ptr(line_ptr).to_str() {
                Ok(s) => line_vec.push(s),
                Err(_) => return ptr::null_mut(),
            }
        }
    }

    match sxs::SxsAlignment::from_lines(&line_vec) {
        Some(aln) => Box::into_raw(Box::new(SxsHandle { aln })),
        None => ptr::null_mut(),
    }
}

/// Free an SXS handle
#[no_mangle]
pub extern "C" fn sxs_free(handle: *mut SxsHandle) {
    if !handle.is_null() {
        unsafe {
            let _ = Box::from_raw(handle);
        }
    }
}

// Field accessors for SXS
#[no_mangle]
pub extern "C" fn sxs_query_sequence_name(handle: *const SxsHandle) -> *mut c_char {
    if handle.is_null() {
        return ptr::null_mut();
    }
    let aln = unsafe { &(*handle).aln };
    match CString::new(aln.query_sequence_name.clone()) {
        Ok(s) => s.into_raw(),
        Err(_) => ptr::null_mut(),
    }
}

#[no_mangle]
pub extern "C" fn sxs_target_sequence_name(handle: *const SxsHandle) -> *mut c_char {
    if handle.is_null() {
        return ptr::null_mut();
    }
    let aln = unsafe { &(*handle).aln };
    match CString::new(aln.target_sequence_name.clone()) {
        Ok(s) => s.into_raw(),
        Err(_) => ptr::null_mut(),
    }
}

#[no_mangle]
pub extern "C" fn sxs_query_start(handle: *const SxsHandle) -> u64 {
    if handle.is_null() {
        return 0;
    }
    unsafe { (*handle).aln.query_start }
}

#[no_mangle]
pub extern "C" fn sxs_query_end(handle: *const SxsHandle) -> u64 {
    if handle.is_null() {
        return 0;
    }
    unsafe { (*handle).aln.query_end }
}

#[no_mangle]
pub extern "C" fn sxs_target_start(handle: *const SxsHandle) -> u64 {
    if handle.is_null() {
        return 0;
    }
    unsafe { (*handle).aln.target_start }
}

#[no_mangle]
pub extern "C" fn sxs_target_end(handle: *const SxsHandle) -> u64 {
    if handle.is_null() {
        return 0;
    }
    unsafe { (*handle).aln.target_end }
}

#[no_mangle]
pub extern "C" fn sxs_num_matches(handle: *const SxsHandle) -> u64 {
    if handle.is_null() {
        return 0;
    }
    unsafe { (*handle).aln.num_matches }
}

#[no_mangle]
pub extern "C" fn sxs_mapping_quality(handle: *const SxsHandle) -> u16 {
    if handle.is_null() {
        return 0;
    }
    unsafe { (*handle).aln.mapping_quality }
}

#[no_mangle]
pub extern "C" fn sxs_cigar(handle: *const SxsHandle) -> *mut CigarHandle {
    if handle.is_null() {
        return ptr::null_mut();
    }
    let aln = unsafe { &(*handle).aln };
    Box::into_raw(Box::new(CigarHandle {
        cigar: aln.cigar.clone(),
    }))
}

#[no_mangle]
pub extern "C" fn sxs_is_good(handle: *const SxsHandle) -> bool {
    if handle.is_null() {
        return false;
    }
    unsafe { (*handle).aln.is_good() }
}

#[no_mangle]
pub extern "C" fn sxs_is_reverse(handle: *const SxsHandle) -> bool {
    if handle.is_null() {
        return false;
    }
    unsafe { (*handle).aln.is_reverse() }
}

// FFI wrappers for compact module

/// Compact nodes by marking boundaries in the graph
///
/// # Arguments
/// * `seqidx_handle` - Handle to the seqindex
/// * `graph_size` - Size of the graph sequence
/// * `node_iitree_handle` - Handle to the node iitree
/// * `path_iitree_handle` - Handle to the path iitree
/// * `seq_id_bv` - Pointer to bitvector array (will be modified)
/// * `seq_id_bv_size` - Size of the bitvector
/// * `num_threads` - Number of threads to use
///
/// # Returns
/// 0 on success, 1 on error
#[no_mangle]
pub extern "C" fn compact_compact_nodes(
    seqidx_handle: *const SeqIndexHandle,
    graph_size: usize,
    node_iitree_handle: *const IITreeHandle,
    path_iitree_handle: *const IITreeHandle,
    seq_id_bv: *mut u64,
    seq_id_bv_size: usize,
    num_threads: usize,
) -> i32 {
    if seqidx_handle.is_null()
        || node_iitree_handle.is_null()
        || path_iitree_handle.is_null()
        || seq_id_bv.is_null()
    {
        eprintln!("[compact] Error: null pointer passed to compact_compact_nodes");
        return 1;
    }

    unsafe {
        let seqidx = Arc::clone(&(*seqidx_handle).seqidx);
        let node_iitree = Arc::clone(&(*node_iitree_handle).iitree);
        let path_iitree = Arc::clone(&(*path_iitree_handle).iitree);

        // Create BitVec from raw pointer
        let bit_count = seq_id_bv_size * 64; // 64 bits per u64
        let slice = std::slice::from_raw_parts_mut(seq_id_bv, seq_id_bv_size);
        let mut bitvec = BitVec::from_slice(slice);
        bitvec.resize(bit_count, false);

        match compact::compact_nodes(
            seqidx,
            graph_size,
            node_iitree,
            path_iitree,
            &mut bitvec,
            num_threads,
        ) {
            Ok(()) => {
                // Copy bitvec back to raw pointer
                let bitvec_slice = bitvec.as_raw_slice();
                std::ptr::copy_nonoverlapping(
                    bitvec_slice.as_ptr(),
                    seq_id_bv,
                    seq_id_bv_size.min(bitvec_slice.len()),
                );
                0
            }
            Err(e) => {
                eprintln!("[compact] Error in compact_nodes: {e}");
                1
            }
        }
    }
}

// FFI wrappers for transclosure module

// TODO: Add helper functions to create handles from C++ objects
// For now, the C++ side needs to manage creation of SeqIndexHandle and IITreeHandle
// directly by wrapping the Rust objects in Arc<> and Arc<Mutex<>>

/// Compute transitive closures for variation graph construction
///
/// # Arguments
/// * `seqidx_handle` - Handle to the seqindex
/// * `aln_iitree_handle` - Handle to the alignment iitree
/// * `seq_v_file` - Path to output sequence file
/// * `node_iitree_handle` - Handle to the node iitree
/// * `path_iitree_handle` - Handle to the path iitree
/// * `repeat_max` - Maximum repeat count
/// * `min_repeat_dist` - Minimum repeat distance
/// * `transclose_batch_size` - Batch size for transitive closure
/// * `show_progress` - Whether to show progress messages
/// * `num_threads` - Number of threads to use
///
/// # Returns
/// The length of the graph sequence, or 0 on error
#[no_mangle]
pub extern "C" fn transclosure_compute(
    seqidx_handle: *const SeqIndexHandle,
    aln_iitree_handle: *const AlnIITreeHandle,
    seq_v_file: *const c_char,
    node_iitree_handle: *const IITreeHandle,
    path_iitree_handle: *const IITreeHandle,
    repeat_max: u64,
    min_repeat_dist: u64,
    transclose_batch_size: u64,
    show_progress: bool,
    num_threads: usize,
) -> usize {
    if seqidx_handle.is_null()
        || aln_iitree_handle.is_null()
        || seq_v_file.is_null()
        || node_iitree_handle.is_null()
        || path_iitree_handle.is_null()
    {
        eprintln!("[transclosure] Error: null pointer passed to transclosure_compute");
        return 0;
    }

    unsafe {
        let seqidx = Arc::clone(&(*seqidx_handle).seqidx);
        let aln_iitree_mutex = Arc::clone(&(*aln_iitree_handle).iitree);
        let node_iitree = Arc::clone(&(*node_iitree_handle).iitree);
        let path_iitree = Arc::clone(&(*path_iitree_handle).iitree);

        // Unwrap the Mutex - alignment tree is read-only during transclosure
        let aln_iitree = match Arc::try_unwrap(aln_iitree_mutex) {
            Ok(mutex) => Arc::new(mutex.into_inner().unwrap()),
            Err(_) => {
                eprintln!("[transclosure] Error: Cannot unwrap aln_iitree Arc (multiple references exist)");
                return 0;
            }
        };

        let seq_v_file_str = match CStr::from_ptr(seq_v_file).to_str() {
            Ok(s) => s,
            Err(e) => {
                eprintln!("[transclosure] Error converting seq_v_file path: {e}");
                return 0;
            }
        };

        match transclosure::compute_transitive_closures(
            seqidx,
            aln_iitree,
            seq_v_file_str,
            node_iitree,
            path_iitree,
            repeat_max,
            min_repeat_dist,
            transclose_batch_size,
            show_progress,
            num_threads,
        ) {
            Ok(length) => length,
            Err(e) => {
                eprintln!("[transclosure] Error in compute_transitive_closures: {e}");
                0
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_add() {
        assert_eq!(seqwish_rust_add(2, 3), 5);
    }
}
