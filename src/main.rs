use clap::Clap;
use rust_htslib::{bam, bam::Read};
use std::fs::File;
use std::io::Write;
use std::path::PathBuf;

#[derive(Clap)]
struct Opts {
    /// Input mapped reads by Bismark aligner (BAM/SAM/CRAM).
    #[clap(parse(from_os_str))]
    input: PathBuf,
    /// Counts of (un)methylated and coverage across cycles (CSV).
    #[clap(parse(from_os_str))]
    mbias: PathBuf,
}

fn main() {
    let opts: Opts = Opts::parse();

    // Create a reader to input file.
    let mut bam = bam::Reader::from_path(&opts.input).unwrap();
    let mut mbias = File::create(opts.mbias).unwrap();

    // Create vectors to store context-aware cytosine methylation status.
    // Context are CpG, CHG and CHH.
    let mut max_len = 0;
    let mut cpg_m = Vec::new();
    let mut cpg_um = Vec::new();
    let mut chg_m = Vec::new();
    let mut chg_um = Vec::new();
    let mut chh_m = Vec::new();
    let mut chh_um = Vec::new();

    // Get statistics about mapped reads.
    let mut total = 0;
    let mut total_valid = 0;

    let mut total_cpg_m = 0;
    let mut total_cpg_um = 0;
    let mut total_chg_m = 0;
    let mut total_chg_um = 0;
    let mut total_chh_m = 0;
    let mut total_chh_um = 0;

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

        // Get Bismark tags about read conversion (XR) and genome conversion (XG)
        let xr = record.aux(b"XR").expect("Missing XR tag").string();
        let xg = record.aux(b"XG").expect("Missing XG tag").string();

        // XR=CT and XG=CT -> original top strand (OT)
        // XR=GA and XG=CT -> complementary to original top strand (CTOT)
        // XR=CT and XG=GA -> original bottom strand (OB)
        // XR=GA and XG=GA -> complementary to original bottom strand (CTOB)
        // TODO: explain why reverse read when CT/GA or GA/CT
        let reverse = (xr == b"CT" && xg == b"GA") || (xr == b"GA" && xg == b"CT");

        // Get CIGAR string as mutable vector.
        let mut cigar = record.cigar().to_vec();

        // Reverse CIGAR string.
        if reverse {
            cigar.reverse();
        }

        // Generate a vector of boolean values according to CIGAR string.
        // Soft-clips (S) and insertions (I) are marked to discard (false),
        // otherwise keep those bases, which are normally matches (M).
        // Each operation is replicated by its length. Example:
        // CIGAR string 3M1D5M3S gives 111011111000.
        let filter: Vec<bool> = cigar
            .iter()
            .flat_map(|c| {
                let keep = c.char() != 'S' || c.char() != 'I';
                vec![keep; c.len() as usize]
            })
            .collect();

        // Get BAM tag `XM` calculated by Bismark aligner.
        // Keep only bases according to vector of boolean values.
        let mut xm: Vec<u8> = record
            .aux(b"XM")
            .unwrap()
            .string()
            .iter()
            .zip(filter)
            .filter(|(_, y)| *y)
            .map(|(x, _)| *x)
            .collect();

        // Reverse read.
        if reverse {
            xm.reverse();
        }

        // Get length of match sequence (without soft-clips or deletions).
        let len = xm.len();

        // If it is greater than action maximun sequence length resise vectors.
        if len > max_len {
            max_len = len;
            cpg_m.resize(max_len, 0);
            cpg_um.resize(max_len, 0);
            chg_m.resize(max_len, 0);
            chg_um.resize(max_len, 0);
            chh_m.resize(max_len, 0);
            chh_um.resize(max_len, 0);
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
                total_chg_m += 1;
            } else if b == b'x' {
                chg_um[i] += 1;
                total_chg_um += 1;
            } else if b == b'Z' {
                cpg_m[i] += 1;
                total_cpg_m += 1;
            } else if b == b'z' {
                cpg_um[i] += 1;
                total_cpg_um += 1;
            } else if b == b'H' {
                chh_m[i] += 1;
                total_chh_m += 1;
            } else if b == b'h' {
                chh_um[i] += 1;
                total_chh_um += 1;
            }
        }
    }

    // Print result tables in CSV format to stdout.
    // Tables are merged as long format and discriminated by `Context` column.
    mbias
        .write_all(b"Context,Cycle,Methylated,Unmethylated,Coverage\n")
        .unwrap();
    for (i, (m, um)) in cpg_m.iter().zip(cpg_um).enumerate() {
        writeln!(&mut mbias, "CpG,{},{},{},{}", i + 1, m, um, m + um).unwrap();
    }
    for (i, (m, um)) in chg_m.iter().zip(chg_um).enumerate() {
        writeln!(&mut mbias, "CHG,{},{},{},{}", i + 1, m, um, m + um).unwrap();
    }
    for (i, (m, um)) in chh_m.iter().zip(chh_um).enumerate() {
        writeln!(&mut mbias, "CHH,{},{},{},{}", i + 1, m, um, m + um).unwrap();
    }

    let percent_valid = total_valid as f32 / total as f32 * 100.0;
    println!("Total:                  {}", total);
    println!(
        "Valid:                  {} ({} %)",
        total_valid, percent_valid
    );

    let total_cpg = total_cpg_m + total_cpg_um;
    let total_chg = total_chg_m + total_chg_um;
    let total_chh = total_chh_m + total_chh_um;
    let total_c = total_cpg + total_chg + total_chh;
    println!("Total C:                {}", total_c);
    println!(
        "Total CpG:              {} ({:.2} %)",
        total_cpg,
        total_cpg as f32 / total_c as f32 * 100.0
    );
    println!(
        "Total CpG Methylated:   {} ({:.2} %)",
        total_cpg_m,
        total_cpg_m as f32 / total_cpg as f32 * 100.0
    );
    println!(
        "Total CpG Unmethylated: {} ({:.2} %)",
        total_cpg_um,
        total_cpg_um as f32 / total_cpg as f32 * 100.0
    );
    println!(
        "Total CHG:              {} ({:.2} %)",
        total_chg,
        total_chg as f32 / total_c as f32 * 100.0
    );
    println!(
        "Total CHG Methylated:   {} ({:.2} %)",
        total_chg_m,
        total_chg_m as f32 / total_chg as f32 * 100.0
    );
    println!(
        "Total CHG Unmethylated: {} ({:.2} %)",
        total_chg_um,
        total_chg_um as f32 / total_chg as f32 * 100.0
    );
    println!(
        "Total CHH:              {} ({:.2} %)",
        total_chh,
        total_chh as f32 / total_c as f32 as f32 * 100.0
    );
    println!(
        "Total CHH Methylated:   {} ({:.2} %)",
        total_chh_m,
        total_chh_m as f32 / total_chh as f32 as f32 * 100.0
    );
    println!(
        "Total CHH Unmethylated: {} ({:.2} %)",
        total_chh_um,
        total_chh_um as f32 / total_chh as f32 as f32 * 100.0
    );
}
