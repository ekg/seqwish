use crate::cigar::{cigar_from_string, cigar_to_string, CigarOp};

/// Represents a single row from a PAF (Pairwise mApping Format) file
#[derive(Debug, Clone, PartialEq)]
pub struct PafRow {
    pub query_sequence_name: String,
    pub query_sequence_length: u64,
    pub query_start: u64,
    pub query_end: u64,
    pub query_target_same_strand: bool,
    pub target_sequence_name: String,
    pub target_sequence_length: u64,
    pub target_start: u64,
    pub target_end: u64,
    pub num_matches: u64,
    pub alignment_block_length: u64,
    pub mapping_quality: u16,
    pub cigar: Vec<CigarOp>,
}

impl PafRow {
    /// Parse a PAF row from a tab/space-delimited line
    pub fn from_line(line: &str) -> Option<Self> {
        let fields: Vec<&str> = line.split([' ', '\t']).collect();

        if fields.len() < 12 {
            return None;
        }

        let query_sequence_name = fields[0].to_string();
        let query_sequence_length = fields[1].parse().ok()?;
        let query_start = fields[2].parse().ok()?;
        let query_end = fields[3].parse().ok()?;
        let query_target_same_strand = fields[4] == "+";
        let target_sequence_name = fields[5].to_string();
        let target_sequence_length = fields[6].parse().ok()?;
        let target_start = fields[7].parse().ok()?;
        let target_end = fields[8].parse().ok()?;
        let num_matches = fields[9].parse().ok()?;
        let alignment_block_length = fields[10].parse().ok()?;
        let mapping_quality = fields[11].parse().ok()?;

        // Find CIGAR in optional fields
        let mut cigar = Vec::new();
        for field in &fields[12..] {
            if let Some(stripped) = field.strip_prefix("cg:Z:") {
                cigar = cigar_from_string(stripped);
                break;
            }
        }

        Some(PafRow {
            query_sequence_name,
            query_sequence_length,
            query_start,
            query_end,
            query_target_same_strand,
            target_sequence_name,
            target_sequence_length,
            target_start,
            target_end,
            num_matches,
            alignment_block_length,
            mapping_quality,
            cigar,
        })
    }

    /// Format PAF row as a tab-delimited string
    #[allow(clippy::inherent_to_string)]
    pub fn to_string(&self) -> String {
        format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\tcg:Z:{}",
            self.query_sequence_name,
            self.query_sequence_length,
            self.query_start,
            self.query_end,
            if self.query_target_same_strand {
                "+"
            } else {
                "-"
            },
            self.target_sequence_name,
            self.target_sequence_length,
            self.target_start,
            self.target_end,
            self.num_matches,
            self.alignment_block_length,
            self.mapping_quality,
            cigar_to_string(&self.cigar)
        )
    }
}

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
        assert_eq!(
            result,
            vec![
                ("file1.paf".to_string(), 0),
                ("file2.paf".to_string(), 0),
                ("file3.paf".to_string(), 0),
            ]
        );
    }

    #[test]
    fn test_parse_paf_spec_mixed() {
        let result = parse_paf_spec("file1.paf:100,file2.paf,file3.paf:250");
        assert_eq!(
            result,
            vec![
                ("file1.paf".to_string(), 100),
                ("file2.paf".to_string(), 0),
                ("file3.paf".to_string(), 250),
            ]
        );
    }

    #[test]
    fn test_parse_paf_spec_empty() {
        let result = parse_paf_spec("");
        assert_eq!(result, Vec::new());
    }

    #[test]
    fn test_parse_paf_spec_whitespace() {
        let result = parse_paf_spec("  file1.paf:100  ,  file2.paf  ");
        assert_eq!(
            result,
            vec![("file1.paf".to_string(), 100), ("file2.paf".to_string(), 0),]
        );
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

    #[test]
    fn test_paf_row_parsing_basic() {
        let line = "query1\t1000\t100\t900\t+\ttarget1\t2000\t200\t1000\t750\t800\t60";
        let row = PafRow::from_line(line).unwrap();

        assert_eq!(row.query_sequence_name, "query1");
        assert_eq!(row.query_sequence_length, 1000);
        assert_eq!(row.query_start, 100);
        assert_eq!(row.query_end, 900);
        assert!(row.query_target_same_strand);
        assert_eq!(row.target_sequence_name, "target1");
        assert_eq!(row.target_sequence_length, 2000);
        assert_eq!(row.target_start, 200);
        assert_eq!(row.target_end, 1000);
        assert_eq!(row.num_matches, 750);
        assert_eq!(row.alignment_block_length, 800);
        assert_eq!(row.mapping_quality, 60);
    }

    #[test]
    fn test_paf_row_parsing_with_cigar() {
        let line =
            "query1\t1000\t100\t900\t-\ttarget1\t2000\t200\t1000\t750\t800\t60\tcg:Z:100M10I50M";
        let row = PafRow::from_line(line).unwrap();

        assert!(!row.query_target_same_strand);
        assert_eq!(row.cigar.len(), 3);
        assert_eq!(row.cigar[0].len, 100);
        assert_eq!(row.cigar[0].op, b'M');
        assert_eq!(row.cigar[1].len, 10);
        assert_eq!(row.cigar[1].op, b'I');
    }

    #[test]
    fn test_paf_row_to_string() {
        let line = "query1\t1000\t100\t900\t+\ttarget1\t2000\t200\t1000\t750\t800\t60\tcg:Z:100M";
        let row = PafRow::from_line(line).unwrap();
        let output = row.to_string();

        assert!(output.contains("query1"));
        assert!(output.contains("target1"));
        assert!(output.contains("cg:Z:100M"));
    }

    #[test]
    fn test_paf_row_parsing_invalid() {
        let line = "query1\t1000\t100"; // Too few fields
        let row = PafRow::from_line(line);
        assert!(row.is_none());
    }
}
