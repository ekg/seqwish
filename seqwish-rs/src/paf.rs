/// Parse a PAF spec string like "file1:100,file2:200,file3"
/// Returns a vector of (filename, weight) pairs
/// If no weight is specified, defaults to 0
pub fn parse_paf_spec(spec: &str) -> Vec<(String, u64)> {
    let mut parsed = Vec::new();

    // Split by comma
    for file_spec in spec.split(',') {
        let file_spec = file_spec.trim();
        if file_spec.is_empty() {
            continue;
        }

        // Split by colon
        let parts: Vec<&str> = file_spec.split(':').collect();

        match parts.len() {
            2 => {
                // filename:weight format
                let filename = parts[0].to_string();
                if let Ok(weight) = parts[1].parse::<u64>() {
                    parsed.push((filename, weight));
                }
            }
            1 => {
                // filename only, weight defaults to 0
                parsed.push((parts[0].to_string(), 0));
            }
            _ => {
                // Invalid format, skip
            }
        }
    }

    parsed
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_paf_spec_single_file() {
        let result = parse_paf_spec("file1.paf");
        assert_eq!(result, vec![("file1.paf".to_string(), 0)]);
    }

    #[test]
    fn test_parse_paf_spec_single_file_with_weight() {
        let result = parse_paf_spec("file1.paf:100");
        assert_eq!(result, vec![("file1.paf".to_string(), 100)]);
    }

    #[test]
    fn test_parse_paf_spec_multiple_files() {
        let result = parse_paf_spec("file1.paf,file2.paf,file3.paf");
        assert_eq!(result, vec![
            ("file1.paf".to_string(), 0),
            ("file2.paf".to_string(), 0),
            ("file3.paf".to_string(), 0),
        ]);
    }

    #[test]
    fn test_parse_paf_spec_mixed() {
        let result = parse_paf_spec("file1.paf:100,file2.paf,file3.paf:250");
        assert_eq!(result, vec![
            ("file1.paf".to_string(), 100),
            ("file2.paf".to_string(), 0),
            ("file3.paf".to_string(), 250),
        ]);
    }

    #[test]
    fn test_parse_paf_spec_empty() {
        let result = parse_paf_spec("");
        assert_eq!(result, Vec::new());
    }

    #[test]
    fn test_parse_paf_spec_whitespace() {
        let result = parse_paf_spec("  file1.paf:100  ,  file2.paf  ");
        assert_eq!(result, vec![
            ("file1.paf".to_string(), 100),
            ("file2.paf".to_string(), 0),
        ]);
    }

    #[test]
    fn test_parse_paf_spec_invalid_weight() {
        let result = parse_paf_spec("file1.paf:abc");
        assert_eq!(result, Vec::new());
    }

    #[test]
    fn test_parse_paf_spec_too_many_colons() {
        let result = parse_paf_spec("file1:100:200");
        assert_eq!(result, Vec::new());
    }
}
