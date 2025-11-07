use seqwish::intervaltree::{AdaptiveTree, IntervalTree};
use std::collections::HashSet;

fn main() -> std::io::Result<()> {
    // Create both tree types
    let mut mem_tree = AdaptiveTree::new_memory();
    let mut disk_tree = AdaptiveTree::new_disk("/tmp/test_tree.iit")?;

    // Add same intervals to both
    mem_tree.open_writer()?;
    disk_tree.open_writer()?;

    // Add some overlapping intervals
    let intervals = vec![
        (10, 20, 1u64),
        (10, 25, 2u64),  // Same start, different end
        (10, 15, 3u64),  // Same start, different end
        (15, 25, 4u64),
        (30, 40, 5u64),
    ];

    for (start, end, val) in &intervals {
        mem_tree.add(*start, *end, *val)?;
        disk_tree.add(*start, *end, *val)?;
    }

    mem_tree.close_writer()?;
    disk_tree.close_writer()?;

    mem_tree.finalize()?;
    disk_tree.finalize()?;

    println!("Added {} intervals", intervals.len());
    println!("\nQuerying range [12, 18):");

    // Query both and collect results
    let mut mem_results = Vec::new();
    let mut disk_results = Vec::new();

    mem_tree.overlap(12, 18, |idx, start, end, val| {
        mem_results.push((idx, start, end, val));
    })?;

    disk_tree.overlap(12, 18, |idx, start, end, val| {
        disk_results.push((idx, start, end, val));
    })?;

    println!("\nMemory results ({}): ", mem_results.len());
    for (idx, start, end, val) in &mem_results {
        println!("  idx={}, [{}, {}), val={}", idx, start, end, val);
    }

    println!("\nDisk results ({}):", disk_results.len());
    for (idx, start, end, val) in &disk_results {
        println!("  idx={}, [{}, {}), val={}", idx, start, end, val);
    }

    // Compare ignoring index (since that might differ due to sort stability)
    let mem_set: HashSet<_> = mem_results.iter().map(|(_, s, e, v)| (s, e, v)).collect();
    let disk_set: HashSet<_> = disk_results.iter().map(|(_, s, e, v)| (s, e, v)).collect();

    if mem_set == disk_set {
        println!("\n✓ Results match (ignoring indices)");
    } else {
        println!("\n✗ Results differ!");
        println!("In memory but not disk: {:?}", mem_set.difference(&disk_set).collect::<Vec<_>>());
        println!("In disk but not memory: {:?}", disk_set.difference(&mem_set).collect::<Vec<_>>());
    }

    Ok(())
}
