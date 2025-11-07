use std::io;
use std::os::unix::io::RawFd;
use std::ptr;

/// Memory-mapped file handle
pub struct MmapHandle {
    pub ptr: *mut libc::c_char,
    pub fd: RawFd,
    pub size: usize,
}

/// Open a file and memory-map it with read-write permissions
///
/// Returns a handle containing the mapped pointer, file descriptor, and file size.
/// The file is mapped with PROT_READ | PROT_WRITE and MAP_SHARED.
/// Memory advice is set to MADV_WILLNEED | MADV_SEQUENTIAL for optimal sequential access.
pub fn mmap_open(filename: &str) -> io::Result<MmapHandle> {
    if filename.is_empty() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "filename is empty",
        ));
    }

    // Open file with read-write permissions
    let c_filename = std::ffi::CString::new(filename)
        .map_err(|_| io::Error::new(io::ErrorKind::InvalidInput, "filename contains null byte"))?;

    let fd = unsafe { libc::open(c_filename.as_ptr(), libc::O_RDWR) };

    if fd == -1 {
        return Err(io::Error::last_os_error());
    }

    // Get file size
    let mut stats: libc::stat = unsafe { std::mem::zeroed() };
    let result = unsafe { libc::fstat(fd, &mut stats) };

    if result == -1 {
        unsafe {
            libc::close(fd);
        }
        return Err(io::Error::last_os_error());
    }

    let fsize = stats.st_size as usize;

    // Memory-map the file
    let ptr = unsafe {
        libc::mmap(
            ptr::null_mut(),
            fsize,
            libc::PROT_READ | libc::PROT_WRITE,
            libc::MAP_SHARED,
            fd,
            0,
        )
    };

    if ptr == libc::MAP_FAILED {
        unsafe {
            libc::close(fd);
        }
        return Err(io::Error::last_os_error());
    }

    // Give memory access hints
    unsafe {
        libc::madvise(ptr, fsize, libc::MADV_WILLNEED | libc::MADV_SEQUENTIAL);
    }

    Ok(MmapHandle {
        ptr: ptr as *mut libc::c_char,
        fd,
        size: fsize,
    })
}

/// Close a memory-mapped file
///
/// Unmaps the memory region and closes the file descriptor.
/// Safe to call multiple times or with already-closed handles.
pub fn mmap_close(handle: &mut MmapHandle) {
    if !handle.ptr.is_null() {
        unsafe {
            libc::munmap(handle.ptr as *mut libc::c_void, handle.size);
        }
        handle.ptr = ptr::null_mut();
    }

    if handle.fd != 0 {
        unsafe {
            libc::close(handle.fd);
        }
        handle.fd = 0;
    }
}

impl Drop for MmapHandle {
    fn drop(&mut self) {
        mmap_close(self);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs;

    #[test]
    fn test_mmap_open_close() {
        // Create a temporary test file
        let test_file = "/tmp/test_mmap_rust.txt";
        let test_data = b"Hello, mmap world!";

        fs::write(test_file, test_data).unwrap();

        // Open and map the file
        let mut handle = mmap_open(test_file).unwrap();
        assert!(!handle.ptr.is_null());
        assert_ne!(handle.fd, 0);
        assert_eq!(handle.size, test_data.len());

        // Read data through mmap
        let mapped_data =
            unsafe { std::slice::from_raw_parts(handle.ptr as *const u8, handle.size) };
        assert_eq!(mapped_data, test_data);

        // Close the mapping
        mmap_close(&mut handle);
        assert!(handle.ptr.is_null());
        assert_eq!(handle.fd, 0);

        // Cleanup
        fs::remove_file(test_file).unwrap();
    }

    #[test]
    fn test_mmap_write_through() {
        let test_file = "/tmp/test_mmap_write_rust.txt";
        let initial_data = b"Initial data here";

        fs::write(test_file, initial_data).unwrap();

        // Open and map the file
        let mut handle = mmap_open(test_file).unwrap();

        // Write through mmap
        unsafe {
            let ptr = handle.ptr as *mut u8;
            ptr::copy_nonoverlapping(b"Modified".as_ptr(), ptr, 8);
        }

        // Sync to ensure write is visible
        unsafe {
            libc::msync(handle.ptr as *mut libc::c_void, handle.size, libc::MS_SYNC);
        }

        mmap_close(&mut handle);

        // Verify the write
        let contents = fs::read(test_file).unwrap();
        assert_eq!(&contents[..8], b"Modified");

        fs::remove_file(test_file).unwrap();
    }

    #[test]
    fn test_mmap_empty_filename() {
        let result = mmap_open("");
        assert!(result.is_err());
    }

    #[test]
    fn test_mmap_nonexistent_file() {
        let result = mmap_open("/tmp/nonexistent_file_12345.txt");
        assert!(result.is_err());
    }

    #[test]
    fn test_mmap_drop_cleanup() {
        let test_file = "/tmp/test_mmap_drop_rust.txt";
        fs::write(test_file, b"Drop test").unwrap();

        {
            let _handle = mmap_open(test_file).unwrap();
            // Handle should be cleaned up when it goes out of scope
        }

        // File should still exist and be accessible
        let contents = fs::read(test_file).unwrap();
        assert_eq!(contents, b"Drop test");

        fs::remove_file(test_file).unwrap();
    }

    #[test]
    fn test_mmap_multiple_close() {
        let test_file = "/tmp/test_mmap_multiclose_rust.txt";
        fs::write(test_file, b"Multi-close test").unwrap();

        let mut handle = mmap_open(test_file).unwrap();

        // First close
        mmap_close(&mut handle);
        assert!(handle.ptr.is_null());

        // Second close should be safe
        mmap_close(&mut handle);
        assert!(handle.ptr.is_null());

        fs::remove_file(test_file).unwrap();
    }
}
