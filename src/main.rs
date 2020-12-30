use rust_htslib::{bam, bam::Read};
use std::env;

fn main() {
    let input = env::args().nth(1).expect("Missing input file.");

    // Get the length of the longest read.
    // Iterate over the entire input file.
    let max_len = bam::Reader::from_path(&input)
        .unwrap()
        .records()
        .map(|record| record.unwrap().seq_len())
        .max()
        .unwrap();

    // Create a reader to input file.
    let mut bam = bam::Reader::from_path(&input).unwrap();

    // Create vectors to store context-aware cytosine methylation status.
    // Context are CpG, CHG and CHH.
    let mut cpg_m = vec![0; max_len];
    let mut cpg_um = vec![0; max_len];
    let mut chg_m = vec![0; max_len];
    let mut chg_um = vec![0; max_len];
    let mut chh_m = vec![0; max_len];
    let mut chh_um = vec![0; max_len];

    // Get statistics about mapped reads.
    let mut total = 0;
    let mut total_valid = 0;
    let mut seq_len = vec![0; max_len];

    // Iterate over the entire input file.
    let mut record = bam::Record::new();
    while let Some(result) = bam.read(&mut record) {
        result.unwrap();
        total += 1;

        // Discard records that are unmapped or duplicated.
        if record.is_unmapped() || record.is_duplicate() {
            continue;
        }
        total_valid += 1;

        // Generate a vector of boolean values according to CIGAR string.
        // Soft-clips (S) and deletions (D) are marked to discard (false),
        // otherwise keep those bases, which are normally matches (M).
        // Each operation is replicated by its length. Example:
        // CIGAR string 3M1D5M3S gives 111011111000.
        let filter: Vec<bool> = record
        .cigar()
            .iter()
            .flat_map(|c| {
                let keep = c.char() != 'S' && c.char() != 'D';
                vec![keep; c.len() as usize]
            })
            .collect();

        // Get BAM tag `XM` calculated by Bismark aligner.
        // Keep only bases according to vector of boolean values.
        let mut xm: Vec<u8> = record.aux(b"XM").unwrap().string()
            .iter()
            .zip(filter)
            .filter(|(_, y)| *y)
            .map(|(x, _)| *x)
            .collect();

        // Count sequence length frequency.
        seq_len[xm.len() - 1] += 1;

        // Get Bismark tags about read conversion (XR) and genome conversion (XG)
        let xr = record.aux(b"XR").expect("Missing XR tag").string();
        let xg = record.aux(b"XG").expect("Missing XG tag").string();

        // XR=CT and XG=CT -> original top strand (OT)
        // XR=GA and XG=CT -> complementary to original top strand (CTOT)
        // XR=CT and XG=GA -> original bottom strand (OB)
        // XR=GA and XG=GA -> complementary to original bottom strand (CTOB)
        // Reverse string.
        if (xr == b"CT" && xg == b"GA") || (xr == b"GA" && xg == b"CT") {
            xm.reverse();
        }

        // Iterate over valid alignment bases (no soft-clips).
        // Count context-aware cytosine methylation states
        //   . for bases not involving cytosines
        //   X for methylated C in CHG context (was protected)
        //   x for not methylated C in CHG context (was converted)
        //   H for methylated C in CHH context (was protected)
        //   h for not methylated C in CHH context (was converted)
        //   Z for methylated C in CpG context (was protected)
        //   z for not methylated C in CpG context (was converted)
        for (i, b) in xm.into_iter().enumerate() {
            if b == b'X' {
                chg_m[i] += 1;
            } else if b == b'x' {
                chg_um[i] += 1;
            } else if b == b'Z' {
                cpg_m[i] += 1;
            } else if b == b'z' {
                cpg_um[i] += 1;
            } else if b == b'H' {
                chh_m[i] += 1;
            } else if b == b'h' {
                chh_um[i] += 1;
            }
        }
    }

    // Print result tables in CSV format to stdout.
    // Tables are merged as long format and discriminated by `Context` column.
    println!("Context,Cycle,Methylated,Unmethylated,Coverage");
    for (i, (m, um)) in cpg_m.iter().zip(cpg_um).enumerate() {
        println!("CpG,{},{},{},{}", i + 1, m, um, m + um);
    }
    for (i, (m, um)) in chg_m.iter().zip(chg_um).enumerate() {
        println!("CHG,{},{},{},{}", i + 1, m, um, m + um);
    }
    for (i, (m, um)) in chh_m.iter().zip(chh_um).enumerate() {
        println!("CHH,{},{},{},{}", i + 1, m, um, m + um);
    }

    let percent_valid = total_valid as f32 / total as f32 * 100.0;
    eprintln!(
        "Total: {}\nValid: {} ({} %)",
        total, total_valid, percent_valid
    );
    eprintln!("Length\tCount");
    for (len, count) in seq_len.into_iter().enumerate() {
        if count > 0 {
            eprintln!("{}\t{}", len + 1, count);
        }
    }
}
