/// CIGAR operation: length + operation character
#[derive(Debug, Clone, PartialEq, Eq)]
#[repr(C)]
pub struct CigarOp {
    pub len: u64,
    pub op: u8,
}

/// Parse CIGAR string into vector of operations
///
/// CIGAR format: alternating runs of digits and operation characters
/// Example: "10M2I5M" -> [(10, 'M'), (2, 'I'), (5, 'M')]
pub fn cigar_from_string(s: &str) -> Vec<CigarOp> {
    // Parse directly over ASCII bytes with an integer accumulator: no String
    // buffer, no per-op reparse. A CIGAR op char follows its digit run, so the
    // op is pushed exactly when its op byte is seen (no trailing flush needed).
    let mut cigar = Vec::with_capacity(s.len() / 2);
    let mut len: u64 = 0;
    let mut have_digit = false;

    for &b in s.as_bytes() {
        if b.is_ascii_digit() {
            len = len * 10 + (b - b'0') as u64;
            have_digit = true;
        } else if have_digit {
            cigar.push(CigarOp { len, op: b });
            len = 0;
            have_digit = false;
        }
    }

    cigar
}

/// Convert CIGAR operations to string
pub fn cigar_to_string(cigar: &[CigarOp]) -> String {
    let mut result = String::new();
    for op in cigar {
        result.push_str(&op.len.to_string());
        result.push(op.op as char);
    }
    result
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_cigar_from_string_simple() {
        let cigar = cigar_from_string("10M");
        assert_eq!(cigar.len(), 1);
        assert_eq!(cigar[0].len, 10);
        assert_eq!(cigar[0].op, b'M');
    }

    #[test]
    fn test_cigar_from_string_multiple_ops() {
        let cigar = cigar_from_string("10M2I5M");
        assert_eq!(cigar.len(), 3);
        assert_eq!(cigar[0], CigarOp { len: 10, op: b'M' });
        assert_eq!(cigar[1], CigarOp { len: 2, op: b'I' });
        assert_eq!(cigar[2], CigarOp { len: 5, op: b'M' });
    }

    #[test]
    fn test_cigar_from_string_complex() {
        let cigar = cigar_from_string("100M1D50M3I25M");
        assert_eq!(cigar.len(), 5);
        assert_eq!(cigar[0], CigarOp { len: 100, op: b'M' });
        assert_eq!(cigar[1], CigarOp { len: 1, op: b'D' });
        assert_eq!(cigar[2], CigarOp { len: 50, op: b'M' });
        assert_eq!(cigar[3], CigarOp { len: 3, op: b'I' });
        assert_eq!(cigar[4], CigarOp { len: 25, op: b'M' });
    }

    #[test]
    fn test_cigar_from_string_various_ops() {
        let cigar = cigar_from_string("5M1I2D3N4S5H6P7X8=");
        assert_eq!(cigar.len(), 9);
        assert_eq!(cigar[0], CigarOp { len: 5, op: b'M' });
        assert_eq!(cigar[1], CigarOp { len: 1, op: b'I' });
        assert_eq!(cigar[2], CigarOp { len: 2, op: b'D' });
        assert_eq!(cigar[3], CigarOp { len: 3, op: b'N' });
        assert_eq!(cigar[4], CigarOp { len: 4, op: b'S' });
        assert_eq!(cigar[5], CigarOp { len: 5, op: b'H' });
        assert_eq!(cigar[6], CigarOp { len: 6, op: b'P' });
        assert_eq!(cigar[7], CigarOp { len: 7, op: b'X' });
        assert_eq!(cigar[8], CigarOp { len: 8, op: b'=' });
    }

    #[test]
    fn test_cigar_from_string_empty() {
        let cigar = cigar_from_string("");
        assert_eq!(cigar.len(), 0);
    }

    #[test]
    fn test_cigar_from_string_large_numbers() {
        let cigar = cigar_from_string("12345678M9876543210D");
        assert_eq!(cigar.len(), 2);
        assert_eq!(
            cigar[0],
            CigarOp {
                len: 12345678,
                op: b'M'
            }
        );
        assert_eq!(
            cigar[1],
            CigarOp {
                len: 9876543210,
                op: b'D'
            }
        );
    }

    #[test]
    fn test_cigar_to_string_simple() {
        let cigar = vec![CigarOp { len: 10, op: b'M' }];
        assert_eq!(cigar_to_string(&cigar), "10M");
    }

    #[test]
    fn test_cigar_to_string_multiple_ops() {
        let cigar = vec![
            CigarOp { len: 10, op: b'M' },
            CigarOp { len: 2, op: b'I' },
            CigarOp { len: 5, op: b'M' },
        ];
        assert_eq!(cigar_to_string(&cigar), "10M2I5M");
    }

    #[test]
    fn test_cigar_to_string_complex() {
        let cigar = vec![
            CigarOp { len: 100, op: b'M' },
            CigarOp { len: 1, op: b'D' },
            CigarOp { len: 50, op: b'M' },
            CigarOp { len: 3, op: b'I' },
            CigarOp { len: 25, op: b'M' },
        ];
        assert_eq!(cigar_to_string(&cigar), "100M1D50M3I25M");
    }

    #[test]
    fn test_cigar_to_string_empty() {
        let cigar = vec![];
        assert_eq!(cigar_to_string(&cigar), "");
    }

    #[test]
    fn test_cigar_roundtrip() {
        let test_cases = vec!["10M", "10M2I5M", "100M1D50M3I25M", "5M1I2D3N4S5H6P7X8=", ""];

        for test in test_cases {
            let cigar = cigar_from_string(test);
            let result = cigar_to_string(&cigar);
            assert_eq!(result, test, "Roundtrip failed for: {}", test);
        }
    }
}
