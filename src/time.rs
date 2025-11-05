use std::time::{SystemTime, UNIX_EPOCH};

/// Get milliseconds since Unix epoch
pub fn time_since_epoch_ms() -> u64 {
    SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .expect("System time before Unix epoch")
        .as_millis() as u64
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_time_since_epoch_ms() {
        let time1 = time_since_epoch_ms();
        // Should be a reasonable timestamp (after year 2020 = 1577836800000 ms)
        assert!(time1 > 1577836800000);

        // Two consecutive calls should be close
        std::thread::sleep(std::time::Duration::from_millis(10));
        let time2 = time_since_epoch_ms();
        assert!(time2 >= time1);
        assert!(time2 - time1 < 1000); // Less than 1 second difference
    }

    #[test]
    fn test_time_monotonically_increasing() {
        let mut times = Vec::new();
        for _ in 0..10 {
            times.push(time_since_epoch_ms());
            std::thread::sleep(std::time::Duration::from_millis(1));
        }

        // Verify times are non-decreasing
        for i in 1..times.len() {
            assert!(times[i] >= times[i-1]);
        }
    }
}
