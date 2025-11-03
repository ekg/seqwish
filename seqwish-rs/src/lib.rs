use std::ffi::{c_char, CStr, CString};
use std::ptr;

pub mod tempfile;

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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_add() {
        assert_eq!(seqwish_rust_add(2, 3), 5);
    }
}
