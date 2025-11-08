use once_cell::sync::Lazy;
use std::collections::HashSet;
use std::ffi::{CStr, CString};
use std::fs;
use std::os::unix::io::FromRawFd;
use std::path::{Path, PathBuf};
use std::sync::Mutex;

/// Global state for temporary file management
struct TempFileState {
    filenames: HashSet<PathBuf>,
    parent_directory: Option<PathBuf>,
    temp_dir: Option<PathBuf>,
    keep_temp: bool,
}

static TEMP_STATE: Lazy<Mutex<TempFileState>> = Lazy::new(|| {
    Mutex::new(TempFileState {
        filenames: HashSet::new(),
        parent_directory: None,
        temp_dir: None,
        keep_temp: false,
    })
});

impl Drop for TempFileState {
    fn drop(&mut self) {
        if !self.keep_temp {
            // Clean up all tracked files
            for filename in &self.filenames {
                let _ = fs::remove_file(filename);
            }

            // Clean up parent directory
            if let Some(ref parent_dir) = self.parent_directory {
                // Remove all remaining files in the directory
                if let Ok(entries) = fs::read_dir(parent_dir) {
                    for entry in entries.flatten() {
                        let _ = fs::remove_file(entry.path());
                    }
                }
                // Remove the directory itself
                let _ = fs::remove_dir(parent_dir);
            }
        }
    }
}

/// Create a temporary file with given base name and suffix
pub fn create(base: &str, suffix: &str) -> Result<PathBuf, std::io::Error> {
    let mut state = TEMP_STATE.lock().unwrap();

    // Create parent directory if needed
    if state.parent_directory.is_none() {
        let temp_dir = get_dir_internal(&state);
        let template = format!("{}/{}", temp_dir.display(), base);

        // Use mkdtemp to create unique directory
        let c_template = CString::new(format!("{template}XXXXXX"))
            .map_err(|_| std::io::Error::new(std::io::ErrorKind::InvalidInput, "invalid path"))?;

        let ptr = unsafe { libc::mkdtemp(c_template.as_ptr() as *mut libc::c_char) };
        if ptr.is_null() {
            return Err(std::io::Error::last_os_error());
        }

        let parent_dir = unsafe { CStr::from_ptr(ptr) }
            .to_str()
            .map_err(|_| std::io::Error::new(std::io::ErrorKind::InvalidInput, "invalid UTF-8"))?
            .to_string();

        state.parent_directory = Some(PathBuf::from(parent_dir));
    }

    // Create the temp file
    let parent = state.parent_directory.as_ref().unwrap();
    let template = format!("{}/XXXXXX{}", parent.display(), suffix);
    let c_template = CString::new(template.clone())
        .map_err(|_| std::io::Error::new(std::io::ErrorKind::InvalidInput, "invalid path"))?;

    let fd = unsafe {
        libc::mkstemps(
            c_template.as_ptr() as *mut libc::c_char,
            suffix.len() as i32,
        )
    };

    if fd == -1 {
        return Err(std::io::Error::last_os_error());
    }

    // Close the file descriptor (we don't keep it open)
    unsafe {
        let file = std::fs::File::from_raw_fd(fd);
        drop(file); // Automatically closes fd
    }

    // Get the actual filename that was created
    let filename = unsafe { CStr::from_ptr(c_template.as_ptr()) }
        .to_str()
        .map_err(|_| std::io::Error::new(std::io::ErrorKind::InvalidInput, "invalid UTF-8"))?
        .to_string();

    let path = PathBuf::from(filename);
    state.filenames.insert(path.clone());

    Ok(path)
}

/// Remove a temporary file
pub fn remove(filename: &Path) {
    let mut state = TEMP_STATE.lock().unwrap();
    let _ = fs::remove_file(filename);
    state.filenames.remove(filename);
}

/// Set the temp directory
pub fn set_dir(new_temp_dir: &str) {
    let mut state = TEMP_STATE.lock().unwrap();
    state.temp_dir = Some(PathBuf::from(new_temp_dir));
}

/// Get the current temp directory
pub fn get_dir() -> PathBuf {
    let state = TEMP_STATE.lock().unwrap();
    get_dir_internal(&state)
}

fn get_dir_internal(state: &TempFileState) -> PathBuf {
    if let Some(ref dir) = state.temp_dir {
        return dir.clone();
    }

    // Use current working directory (matches C++ behavior)
    std::env::current_dir().unwrap_or_else(|_| PathBuf::from("."))
}

/// Set whether to keep temp files on exit
pub fn set_keep_temp(setting: bool) {
    let mut state = TEMP_STATE.lock().unwrap();
    state.keep_temp = setting;
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs;

    #[test]
    fn test_create_and_remove() {
        let path = create("test", ".tmp").expect("Failed to create temp file");
        assert!(path.exists(), "Temp file should exist");
        assert!(
            path.to_string_lossy().ends_with(".tmp"),
            "Should have correct suffix"
        );

        remove(&path);
        assert!(!path.exists(), "Temp file should be removed");
    }

    #[test]
    fn test_set_get_dir() {
        let original = get_dir();

        set_dir("/tmp");
        assert_eq!(get_dir(), PathBuf::from("/tmp"));

        // Restore original
        set_dir(&original.to_string_lossy());
    }

    #[test]
    fn test_multiple_files() {
        let file1 = create("test", ".txt").expect("Failed to create file1");
        let file2 = create("test", ".txt").expect("Failed to create file2");

        assert!(file1.exists());
        assert!(file2.exists());
        assert_ne!(file1, file2, "Files should have unique names");

        remove(&file1);
        remove(&file2);
    }

    #[test]
    fn test_file_is_writable() {
        let path = create("test", ".data").expect("Failed to create temp file");

        // Try to write to it
        fs::write(&path, b"hello world").expect("Should be able to write");
        let content = fs::read(&path).expect("Should be able to read");
        assert_eq!(content, b"hello world");

        remove(&path);
    }
}
