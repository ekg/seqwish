use std::ffi::CString;

/// Check if a file exists using stat
pub fn file_exists(filename: &str) -> bool {
    let c_filename = match CString::new(filename) {
        Ok(s) => s,
        Err(_) => return false,
    };

    let mut stat_buf: libc::stat = unsafe { std::mem::zeroed() };
    let result = unsafe { libc::stat(c_filename.as_ptr(), &mut stat_buf) };

    result == 0
}

/// Parse a number with optional suffix (k, m, g for thousands, millions, billions)
/// Returns default_value if parsing fails
pub fn handy_parameter(value: &str, default_value: f64) -> f64 {
    fn is_a_number(s: &str) -> bool {
        if s.is_empty() {
            return false;
        }

        let dot_count = s.chars().filter(|&c| c == '.').count();
        if dot_count >= 2 {
            return false;
        }

        s.chars().all(|c| c.is_ascii_digit() || c == '.')
    }

    let mut str_len = value.len();
    let mut exp: u32 = 0;

    if str_len > 0 {
        match value.chars().last() {
            Some('k') | Some('K') => {
                exp = 3;
                str_len -= 1;
            }
            Some('m') | Some('M') => {
                exp = 6;
                str_len -= 1;
            }
            Some('g') | Some('G') => {
                exp = 9;
                str_len -= 1;
            }
            _ => {}
        }
    }

    let tmp = &value[..str_len];

    if is_a_number(tmp) {
        match tmp.parse::<f64>() {
            Ok(num) => num * 10_f64.powi(exp as i32),
            Err(_) => default_value,
        }
    } else {
        default_value
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs;

    #[test]
    fn test_file_exists_true() {
        let test_file = "/tmp/test_file_exists_rust.txt";
        fs::write(test_file, b"test").unwrap();

        assert!(file_exists(test_file));

        fs::remove_file(test_file).unwrap();
    }

    #[test]
    fn test_file_exists_false() {
        assert!(!file_exists("/tmp/nonexistent_file_12345678.txt"));
    }

    #[test]
    fn test_file_exists_empty_string() {
        assert!(!file_exists(""));
    }

    #[test]
    fn test_handy_parameter_plain_number() {
        assert_eq!(handy_parameter("42", 0.0), 42.0);
        assert_eq!(handy_parameter("3.15", 0.0), 3.15);
    }

    #[test]
    fn test_handy_parameter_with_k() {
        assert_eq!(handy_parameter("10k", 0.0), 10_000.0);
        assert_eq!(handy_parameter("10K", 0.0), 10_000.0);
        assert_eq!(handy_parameter("1.5k", 0.0), 1_500.0);
    }

    #[test]
    fn test_handy_parameter_with_m() {
        assert_eq!(handy_parameter("10m", 0.0), 10_000_000.0);
        assert_eq!(handy_parameter("10M", 0.0), 10_000_000.0);
        assert_eq!(handy_parameter("2.5m", 0.0), 2_500_000.0);
    }

    #[test]
    fn test_handy_parameter_with_g() {
        assert_eq!(handy_parameter("10g", 0.0), 10_000_000_000.0);
        assert_eq!(handy_parameter("10G", 0.0), 10_000_000_000.0);
        assert_eq!(handy_parameter("1.2g", 0.0), 1_200_000_000.0);
    }

    #[test]
    fn test_handy_parameter_invalid() {
        assert_eq!(handy_parameter("abc", 99.0), 99.0);
        assert_eq!(handy_parameter("", 99.0), 99.0);
        assert_eq!(handy_parameter("1.2.3", 99.0), 99.0);
    }

    #[test]
    fn test_handy_parameter_no_suffix() {
        assert_eq!(handy_parameter("100", 50.0), 100.0);
    }

    #[test]
    fn test_handy_parameter_zero() {
        assert_eq!(handy_parameter("0", 99.0), 0.0);
        assert_eq!(handy_parameter("0k", 99.0), 0.0);
    }
}
