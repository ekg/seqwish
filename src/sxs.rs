use crate::cigar::{cigar_from_string, cigar_to_string, CigarOp};

/// Represents an SXS (Seqwish eXtended format) alignment
#[derive(Debug, Clone, PartialEq)]
pub struct SxsAlignment {
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
    pub mapping_quality: u16,
    pub cigar: Vec<CigarOp>,
}

impl SxsAlignment {
    /// Create a new empty SXS alignment
    pub fn new() -> Self {
        SxsAlignment {
            query_sequence_name: String::new(),
            query_sequence_length: 0,
            query_start: 0,
            query_end: 0,
            query_target_same_strand: true,
            target_sequence_name: String::new(),
            target_sequence_length: 0,
            target_start: 0,
            target_end: 0,
            num_matches: 0,
            mapping_quality: 0,
            cigar: Vec::new(),
        }
    }

    /// Parse an SXS alignment from lines
    /// Lines should be in the format:
    /// A <target_name> <query_name>
    /// I <target_start> <target_end> <query_start> <query_end>
    /// M <num_matches>
    /// C <cigar_string>
    /// Q <mapping_quality>
    pub fn from_lines(lines: &[&str]) -> Option<Self> {
        let mut aln = SxsAlignment::new();

        for line in lines {
            let line = line.trim();
            if line.is_empty() {
                continue;
            }

            let fields: Vec<&str> = line.split(|c| c == ' ' || c == '\t').collect();
            if fields.is_empty() {
                continue;
            }

            match fields[0] {
                "A" => {
                    if fields.len() >= 3 {
                        aln.target_sequence_name = fields[1].to_string();
                        aln.query_sequence_name = fields[2].to_string();
                    }
                }
                "I" => {
                    if fields.len() >= 5 {
                        aln.target_start = fields[1].parse().ok()?;
                        aln.target_end = fields[2].parse().ok()?;
                        aln.query_start = fields[3].parse().ok()?;
                        aln.query_end = fields[4].parse().ok()?;
                    }
                }
                "M" => {
                    if fields.len() >= 2 {
                        aln.num_matches = fields[1].parse().ok()?;
                    }
                }
                "C" => {
                    if fields.len() >= 2 {
                        aln.cigar = cigar_from_string(fields[1]);
                    }
                }
                "Q" => {
                    if fields.len() >= 2 {
                        aln.mapping_quality = fields[1].parse().ok()?;
                    }
                }
                "T" => {
                    // unhandled, skip
                }
                _ => {
                    // unknown line type, skip
                }
            }
        }

        // Check if we have at least sequence names
        if aln.query_sequence_name.is_empty() {
            return None;
        }

        Some(aln)
    }

    /// Check if alignment is valid (has sequence names)
    pub fn is_good(&self) -> bool {
        !self.query_sequence_name.is_empty()
    }

    /// Check if query is reverse complement
    pub fn is_reverse(&self) -> bool {
        self.query_start > self.query_end
    }

    /// Format as SXS string
    pub fn to_string(&self) -> String {
        format!(
            "A\t{}\t{}\nI\t{}\t{}\t{}\t{}\nM\t{}\nC\t{}\nQ\t{}",
            self.target_sequence_name,
            self.query_sequence_name,
            self.target_start,
            self.target_end,
            self.query_start,
            self.query_end,
            self.num_matches,
            cigar_to_string(&self.cigar),
            self.mapping_quality
        )
    }
}

impl Default for SxsAlignment {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sxs_parsing_basic() {
        let lines = vec![
            "A\ttarget1\tquery1",
            "I\t100\t200\t50\t150",
            "M\t90",
            "C\t100M",
            "Q\t60",
        ];

        let aln = SxsAlignment::from_lines(&lines).unwrap();

        assert_eq!(aln.target_sequence_name, "target1");
        assert_eq!(aln.query_sequence_name, "query1");
        assert_eq!(aln.target_start, 100);
        assert_eq!(aln.target_end, 200);
        assert_eq!(aln.query_start, 50);
        assert_eq!(aln.query_end, 150);
        assert_eq!(aln.num_matches, 90);
        assert_eq!(aln.mapping_quality, 60);
        assert_eq!(aln.cigar.len(), 1);
        assert_eq!(aln.cigar[0].len, 100);
        assert_eq!(aln.cigar[0].op, b'M');
    }

    #[test]
    fn test_sxs_is_good() {
        let mut aln = SxsAlignment::new();
        assert!(!aln.is_good());

        aln.query_sequence_name = "query1".to_string();
        assert!(aln.is_good());
    }

    #[test]
    fn test_sxs_is_reverse() {
        let mut aln = SxsAlignment::new();
        aln.query_start = 100;
        aln.query_end = 50;
        assert!(aln.is_reverse());

        aln.query_start = 50;
        aln.query_end = 100;
        assert!(!aln.is_reverse());
    }

    #[test]
    fn test_sxs_to_string() {
        let lines = vec![
            "A\ttarget1\tquery1",
            "I\t100\t200\t50\t150",
            "M\t90",
            "C\t100M",
            "Q\t60",
        ];

        let aln = SxsAlignment::from_lines(&lines).unwrap();
        let output = aln.to_string();

        assert!(output.contains("target1"));
        assert!(output.contains("query1"));
        assert!(output.contains("100M"));
    }

    #[test]
    fn test_sxs_parsing_with_t_line() {
        let lines = vec![
            "A\ttarget1\tquery1",
            "I\t100\t200\t50\t150",
            "T\tignored_data",
            "M\t90",
            "C\t50M10I40M",
            "Q\t60",
        ];

        let aln = SxsAlignment::from_lines(&lines).unwrap();

        assert_eq!(aln.target_sequence_name, "target1");
        assert_eq!(aln.cigar.len(), 3);
    }

    #[test]
    fn test_sxs_parsing_invalid() {
        let lines = vec!["I\t100\t200\t50\t150", "M\t90"]; // No 'A' line

        let aln = SxsAlignment::from_lines(&lines);
        assert!(aln.is_none());
    }
}
