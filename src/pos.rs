/// Position type - encodes offset and orientation in a single u64
/// Bottom bit: orientation (0 = forward, 1 = reverse)
/// Upper 63 bits: offset in sequence
pub type PosT = u64;

/// Position with alignment length
#[repr(C)]
#[derive(Debug, Copy, Clone, PartialEq, Eq)]
pub struct AlnPosT {
    pub pos: PosT,
    pub aln_length: u64,
}

impl PartialOrd for AlnPosT {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for AlnPosT {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        // Match C++ behavior: a < b if both pos and aln_length are less
        match self.pos.cmp(&other.pos) {
            std::cmp::Ordering::Less => {
                if self.aln_length < other.aln_length {
                    std::cmp::Ordering::Less
                } else {
                    std::cmp::Ordering::Equal
                }
            }
            std::cmp::Ordering::Equal => self.aln_length.cmp(&other.aln_length),
            std::cmp::Ordering::Greater => std::cmp::Ordering::Greater,
        }
    }
}

/// Create a position from offset and orientation
#[inline]
pub fn make_pos_t(offset: u64, is_rev: bool) -> PosT {
    // Top bit is reserved for is_rev flag
    // The rest is our offset in the input sequence vector
    let rev_mask: u64 = 1;
    let mut pos = offset << 1;
    // Conditional set/clear without branching
    pos = (pos & !rev_mask) | ((is_rev as u64).wrapping_neg() & rev_mask);
    pos
}

/// Extract offset from position
#[inline]
pub fn offset(pos: PosT) -> u64 {
    pos >> 1
}

/// Check if position is reverse orientation
#[inline]
pub fn is_rev(pos: PosT) -> bool {
    (pos & 1) != 0
}

/// Increment position (moves forward in sequence orientation)
#[inline]
pub fn incr_pos(pos: &mut PosT) {
    if is_rev(*pos) {
        *pos -= 2;
    } else {
        *pos += 2;
    }
}

/// Increment position by N bases
#[inline]
pub fn incr_pos_by(pos: &mut PosT, by: usize) {
    if is_rev(*pos) {
        *pos -= 2 * by as u64;
    } else {
        *pos += 2 * by as u64;
    }
}

/// Decrement position (moves backward in sequence orientation)
/// Note: Uses wrapping subtraction to handle underflow (when pos < 2)
/// This matches C++ behavior where underflow wraps to a large number
#[inline]
pub fn decr_pos(pos: &mut PosT) {
    if !is_rev(*pos) {
        *pos = pos.wrapping_sub(2);
    } else {
        *pos += 2;
    }
}

/// Decrement position by N bases
/// Note: Uses wrapping subtraction to handle underflow
/// This matches C++ behavior where underflow wraps to a large number
#[inline]
pub fn decr_pos_by(pos: &mut PosT, by: usize) {
    if !is_rev(*pos) {
        *pos = pos.wrapping_sub(2 * by as u64);
    } else {
        *pos += 2 * by as u64;
    }
}

/// Reverse the orientation of a position
#[inline]
pub fn rev_pos_t(pos: PosT) -> PosT {
    make_pos_t(offset(pos), !is_rev(pos))
}

/// Convert position to string representation (e.g., "42+" or "123-")
pub fn pos_to_string(pos: PosT) -> String {
    format!("{}{}", offset(pos), if is_rev(pos) { "-" } else { "+" })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_make_pos_t() {
        // Forward position at offset 0
        let pos = make_pos_t(0, false);
        assert_eq!(pos, 0);
        assert_eq!(offset(pos), 0);
        assert!(!is_rev(pos));

        // Reverse position at offset 0
        let pos = make_pos_t(0, true);
        assert_eq!(pos, 1);
        assert_eq!(offset(pos), 0);
        assert!(is_rev(pos));

        // Forward position at offset 100
        let pos = make_pos_t(100, false);
        assert_eq!(pos, 200);
        assert_eq!(offset(pos), 100);
        assert!(!is_rev(pos));

        // Reverse position at offset 100
        let pos = make_pos_t(100, true);
        assert_eq!(pos, 201);
        assert_eq!(offset(pos), 100);
        assert!(is_rev(pos));
    }

    #[test]
    fn test_offset_extraction() {
        assert_eq!(offset(make_pos_t(0, false)), 0);
        assert_eq!(offset(make_pos_t(0, true)), 0);
        assert_eq!(offset(make_pos_t(1000, false)), 1000);
        assert_eq!(offset(make_pos_t(1000, true)), 1000);
        assert_eq!(offset(make_pos_t(u64::MAX >> 1, false)), u64::MAX >> 1);
    }

    #[test]
    fn test_is_rev() {
        assert!(!is_rev(make_pos_t(0, false)));
        assert!(is_rev(make_pos_t(0, true)));
        assert!(!is_rev(make_pos_t(12345, false)));
        assert!(is_rev(make_pos_t(12345, true)));
    }

    #[test]
    fn test_incr_pos() {
        // Forward: increment increases encoded value
        let mut pos = make_pos_t(10, false);
        incr_pos(&mut pos);
        assert_eq!(offset(pos), 11);
        assert!(!is_rev(pos));

        // Reverse: increment decreases encoded value
        let mut pos = make_pos_t(10, true);
        incr_pos(&mut pos);
        assert_eq!(offset(pos), 9);
        assert!(is_rev(pos));
    }

    #[test]
    fn test_incr_pos_by() {
        let mut pos = make_pos_t(100, false);
        incr_pos_by(&mut pos, 5);
        assert_eq!(offset(pos), 105);

        let mut pos = make_pos_t(100, true);
        incr_pos_by(&mut pos, 5);
        assert_eq!(offset(pos), 95);
    }

    #[test]
    fn test_decr_pos() {
        // Forward: decrement decreases encoded value
        let mut pos = make_pos_t(10, false);
        decr_pos(&mut pos);
        assert_eq!(offset(pos), 9);
        assert!(!is_rev(pos));

        // Reverse: decrement increases encoded value
        let mut pos = make_pos_t(10, true);
        decr_pos(&mut pos);
        assert_eq!(offset(pos), 11);
        assert!(is_rev(pos));
    }

    #[test]
    fn test_decr_pos_by() {
        let mut pos = make_pos_t(100, false);
        decr_pos_by(&mut pos, 5);
        assert_eq!(offset(pos), 95);

        let mut pos = make_pos_t(100, true);
        decr_pos_by(&mut pos, 5);
        assert_eq!(offset(pos), 105);
    }

    #[test]
    fn test_rev_pos_t() {
        let pos_fwd = make_pos_t(42, false);
        let pos_rev = rev_pos_t(pos_fwd);
        assert_eq!(offset(pos_rev), 42);
        assert!(is_rev(pos_rev));

        let pos_fwd_again = rev_pos_t(pos_rev);
        assert_eq!(pos_fwd, pos_fwd_again);
    }

    #[test]
    fn test_pos_to_string() {
        assert_eq!(pos_to_string(make_pos_t(0, false)), "0+");
        assert_eq!(pos_to_string(make_pos_t(0, true)), "0-");
        assert_eq!(pos_to_string(make_pos_t(42, false)), "42+");
        assert_eq!(pos_to_string(make_pos_t(123, true)), "123-");
    }

    #[test]
    fn test_aln_pos_t_equality() {
        let a = AlnPosT {
            pos: make_pos_t(10, false),
            aln_length: 100,
        };
        let b = AlnPosT {
            pos: make_pos_t(10, false),
            aln_length: 100,
        };
        let c = AlnPosT {
            pos: make_pos_t(10, false),
            aln_length: 200,
        };

        assert_eq!(a, b);
        assert_ne!(a, c);
    }

    #[test]
    fn test_aln_pos_t_ordering() {
        let a = AlnPosT {
            pos: make_pos_t(5, false),
            aln_length: 50,
        };
        let b = AlnPosT {
            pos: make_pos_t(10, false),
            aln_length: 100,
        };
        let c = AlnPosT {
            pos: make_pos_t(10, false),
            aln_length: 50,
        };

        // a has both smaller pos and aln_length than b
        assert!(a < b);

        // c has smaller aln_length but same pos as b
        assert!(c < b);
    }

    #[test]
    fn test_roundtrip() {
        for offset_val in [0u64, 1, 100, 1000, 1_000_000, u64::MAX >> 1] {
            for &is_rev_val in &[false, true] {
                let pos = make_pos_t(offset_val, is_rev_val);
                assert_eq!(offset(pos), offset_val);
                assert_eq!(is_rev(pos), is_rev_val);
            }
        }
    }
}
