use std::collections::HashMap;
use std::ffi::CString;
use std::os::raw::c_char;
use once_cell::sync::Lazy;

// Git version injected at build time - fallback to "unknown" if not set
const GIT_VERSION: &str = match option_env!("SEQWISH_GIT_VERSION") {
    Some(v) => v,
    None => "unknown",
};

static CODENAMES: Lazy<HashMap<&'static str, &'static str>> = Lazy::new(|| {
    let mut m = HashMap::new();
        m.insert("v0.1", "initial release");
        m.insert("v0.2", "parallel union find / transclosure");
        m.insert("v0.2.1", "unstable point release");
        m.insert("v0.4", "improving stability");
        m.insert("v0.4.1", "complete correction of range compression in union find");
        m.insert("v0.5", "repeat limitation");
        m.insert("0.6", "std::thread our way through things");
        m.insert("v0.6.1", "simplify paths in GFA");
        m.insert("v0.7", "memory safety and minor improvements");
        m.insert("v0.7.1", "seqwish 0.7.1 - Sicurezza");
        m.insert("v0.7.2", "seqwish 0.7.2 - Diffidenza");
        m.insert("v0.7.3", "seqwish 0.7.3 - Veloce");
        m.insert("v0.7.4", "seqwish 0.7.4 - Rifinitura");
        m.insert("v0.7.5", "seqwish 0.7.5 - Pasticcione");
        m.insert("v0.7.6", "seqwish 0.7.6 - Temporaneo");
        m
});

pub fn get_version() -> &'static str {
    GIT_VERSION
}

pub fn get_release() -> String {
    match GIT_VERSION.find('-') {
        None => GIT_VERSION.to_string(),
        Some(dash_pos) => GIT_VERSION[..dash_pos].to_string(),
    }
}

pub fn get_codename() -> Option<&'static str> {
    let release = get_release();
    CODENAMES.get(release.as_str()).copied()
}

pub fn get_short() -> String {
    match get_codename() {
        Some(codename) => format!("{} \"{}\"", GIT_VERSION, codename),
        None => GIT_VERSION.to_string(),
    }
}

// FFI exports for C++ compatibility
#[no_mangle]
pub extern "C" fn version_get_version() -> *mut c_char {
    CString::new(get_version()).unwrap().into_raw()
}

#[no_mangle]
pub extern "C" fn version_get_release() -> *mut c_char {
    CString::new(get_release()).unwrap().into_raw()
}

#[no_mangle]
pub extern "C" fn version_get_codename() -> *mut c_char {
    match get_codename() {
        Some(codename) => CString::new(codename).unwrap().into_raw(),
        None => std::ptr::null_mut(),
    }
}

#[no_mangle]
pub extern "C" fn version_get_short() -> *mut c_char {
    CString::new(get_short()).unwrap().into_raw()
}

#[no_mangle]
pub extern "C" fn version_free_string(s: *mut c_char) {
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
    fn test_get_version() {
        let version = get_version();
        assert!(!version.is_empty());
    }

    #[test]
    fn test_get_release() {
        // Should return the part before the dash, or the whole string if no dash
        let release = get_release();
        assert!(!release.is_empty());
    }

    #[test]
    fn test_get_short() {
        let short = get_short();
        assert!(!short.is_empty());
        assert!(short.contains(GIT_VERSION));
    }

    #[test]
    fn test_codenames_exist() {
        assert!(CODENAMES.contains_key("v0.7.6"));
        assert_eq!(CODENAMES.get("v0.7.6"), Some(&"seqwish 0.7.6 - Temporaneo"));
    }
}
