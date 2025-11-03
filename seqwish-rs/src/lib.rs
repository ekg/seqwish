use std::ffi::{c_char, CStr, CString};
use std::ptr;

pub mod tempfile;
pub mod pos;
pub mod dna;
pub mod cigar;
pub mod mmap;

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

/// Parse CIGAR string and return handle to CIGAR vector
/// Returns NULL on error. Must be freed with cigar_free.
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
pub extern "C" fn cigar_get_op(handle: *const CigarHandle, index: usize, len_out: *mut u64, op_out: *mut u8) -> bool {
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
pub extern "C" fn cigar_free(handle: *mut CigarHandle) {
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

    let mut handle = mmap::MmapHandle {
        ptr: buf,
        fd,
        size,
    };

    mmap::mmap_close(&mut handle);
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_add() {
        assert_eq!(seqwish_rust_add(2, 3), 5);
    }
}
