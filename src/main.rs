use clap::Clap;
use rust_htslib::bam::record::Cigar;
use rust_htslib::{bam, bam::Read};
use std::collections::btree_map::BTreeMap;
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
    /// Counts of (un)methylated and coverage across genome (cov).
    #[clap(parse(from_os_str))]
    bismark_cov: PathBuf,
    /// Counts of (un)methylated and coverage across genome (bedGraph).
    #[clap(parse(from_os_str))]
    bed_graph: PathBuf,
    /// Starting vector capacity to store per sequencing cycles data.
    #[clap(long, default_value = "200")]
    start_capacity: usize,
}

fn main() {
    let opts: Opts = Opts::parse();

    // Create a reader to input file.
    let mut bam = bam::Reader::from_path(&opts.input).unwrap();
    let mut mbias = File::create(opts.mbias).unwrap();
    let mut bismark_cov = File::create(opts.bismark_cov).unwrap();
    let mut bed_graph = File::create(opts.bed_graph).unwrap();

    let header = bam::Header::from_template(bam.header()).to_hashmap();

    // Create vectors to store context-aware cytosine methylation status.
    // Context are CpG, CHG and CHH.
    let mut max_len = 0;
    let mut read_cpg_m = Vec::with_capacity(opts.start_capacity);
    let mut read_cpg_um = Vec::with_capacity(opts.start_capacity);
    let mut read_chg_m = Vec::with_capacity(opts.start_capacity);
    let mut read_chg_um = Vec::with_capacity(opts.start_capacity);
    let mut read_chh_m = Vec::with_capacity(opts.start_capacity);
    let mut read_chh_um = Vec::with_capacity(opts.start_capacity);

    // Get statistics about mapped reads.
    let mut total = 0;
    let mut total_valid = 0;

    let mut total_cpg_m = 0;
    let mut total_cpg_um = 0;
    let mut total_chg_m = 0;
    let mut total_chg_um = 0;
    let mut total_chh_m = 0;
    let mut total_chh_um = 0;

    // Methylation calls across genome.
    // Methylated, unmethylated and coverage.
    let mut chr_cpg: BTreeMap<String, BTreeMap<i64, (u32, u32, u32)>> = BTreeMap::new();
    let mut chr_chg: BTreeMap<String, BTreeMap<i64, (u32, u32, u32)>> = BTreeMap::new();
    let mut chr_chh: BTreeMap<String, BTreeMap<i64, (u32, u32, u32)>> = BTreeMap::new();

    // Iterate over SQ header tags to initialize chromosomes.
    let contigs = header.get("SQ").unwrap();
    for contig in contigs {
        let chr = contig.get("SN").unwrap();
        chr_cpg.insert(chr.to_string(), BTreeMap::new());
        chr_chg.insert(chr.to_string(), BTreeMap::new());
        chr_chh.insert(chr.to_string(), BTreeMap::new());
    }

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
        let cigar = record.cigar();

        // Get BAM tag `XM` calculated by Bismark aligner.
        // Keep only bases according to vector of boolean values.
        let xm = record.aux(b"XM").unwrap().string().to_owned();

        // Get length of match sequence (without soft-clips or insertions).
        let len = cigar
            .iter()
            .filter(|&c| c.char() != 'S' || c.char() != 'I')
            .map(|&c| c.len() as usize)
            .sum();

        // If it is greater than action maximun sequence length resise vectors.
        if len > max_len {
            max_len = len;
            read_cpg_m.resize(max_len, 0);
            read_cpg_um.resize(max_len, 0);
            read_chg_m.resize(max_len, 0);
            read_chg_um.resize(max_len, 0);
            read_chh_m.resize(max_len, 0);
            read_chh_um.resize(max_len, 0);
        }

        let chr = bam.header().tid2name(record.tid() as u32).to_owned();
        let chr = String::from_utf8(chr).unwrap();
        let cpg = chr_cpg.get_mut(&chr).unwrap();
        let chg = chr_chg.get_mut(&chr).unwrap();
        let chh = chr_chh.get_mut(&chr).unwrap();

        // Iterate over valid alignment bases (no soft-clips).
        // Count context-aware cytosine methylation states
        //   . for bases not involving cytosines
        //   X for methylated C in CHG context (was protected)
        //   x for not methylated C in CHG context (was converted)
        //   H for methylated C in CHH context (was protected)
        //   h for not methylated C in CHH context (was converted)
        //   Z for methylated C in CpG context (was protected)
        //   z for not methylated C in CpG context (was converted)
        let mut xm_idx: usize = 0;
        let mut read_idx: usize = if !reverse { 0 } else { len };
        let mut pos = record.pos();
        for c in cigar.iter() {
            match c {
                Cigar::Match(len) => {
                    while xm_idx < *len as usize {
                        if reverse {
                            read_idx -= 1;
                        }
                        let b = xm[xm_idx];
                        if b == b'X' || b == b'x' {
                            let (m, um, cov) = chg.entry(pos).or_insert((0, 0, 0));
                            if b == b'X' {
                                read_chg_m[read_idx] += 1;
                                total_chg_m += 1;
                                *m += 1;
                            } else {
                                read_chg_um[read_idx] += 1;
                                total_chg_um += 1;
                                *um += 1;
                            }
                            *cov += 1;
                        } else if b == b'Z' || b == b'z' {
                            let (m, um, cov) = cpg.entry(pos).or_insert((0, 0, 0));
                            if b == b'Z' {
                                read_cpg_m[read_idx] += 1;
                                total_cpg_m += 1;
                                *m += 1;
                            } else {
                                read_cpg_um[read_idx] += 1;
                                total_cpg_um += 1;
                                *um += 1;
                            }
                            *cov += 1;
                        } else if b == b'H' || b == b'h' {
                            let (m, um, cov) = chh.entry(pos).or_insert((0, 0, 0));
                            if b == b'H' {
                                read_chh_m[read_idx] += 1;
                                total_chh_m += 1;
                                *m += 1;
                            } else {
                                read_chh_um[read_idx] += 1;
                                total_chh_um += 1;
                                *um += 1;
                            }
                            *cov += 1;
                        }
                        pos += 1;
                        xm_idx += 1;
                        if !reverse {
                            read_idx += 1;
                        }
                    }
                }
                Cigar::Ins(len) => {
                    xm_idx += *len as usize;
                    if !reverse {
                        read_idx += *len as usize;
                    } else {
                        read_idx -= *len as usize;
                    }
                }
                Cigar::Del(len) => {
                    pos += *len as i64;
                }
                Cigar::SoftClip(len) => {
                    xm_idx += *len as usize;
                }
                _ => {
                    panic!("Invalid CIGAR operation: {}", c.char());
                }
            }
        }
    }

    // Print result tables in CSV format to stdout.
    // Tables are merged as long format and discriminated by `Context` column.
    mbias
        .write_all(b"Context,Cycle,Methylated,Unmethylated,Coverage\n")
        .unwrap();
    for (i, (m, um)) in read_cpg_m.iter().zip(read_cpg_um).enumerate() {
        writeln!(&mut mbias, "CpG,{},{},{},{}", i + 1, m, um, m + um).unwrap();
    }
    for (i, (m, um)) in read_chg_m.iter().zip(read_chg_um).enumerate() {
        writeln!(&mut mbias, "CHG,{},{},{},{}", i + 1, m, um, m + um).unwrap();
    }
    for (i, (m, um)) in read_chh_m.iter().zip(read_chh_um).enumerate() {
        writeln!(&mut mbias, "CHH,{},{},{},{}", i + 1, m, um, m + um).unwrap();
    }

    bed_graph.write_all(b"track type=bedGraph\n").unwrap();
    for (chr, xs) in chr_cpg {
        for (pos, (m, u, cov)) in xs {
            let perc = m as f64 / cov as f64 * 100.0;
            let pos1 = pos + 1;
            writeln!(
                &mut bismark_cov,
                "{}\t{}\t{}\t{}\t{}\t{}",
                chr, pos1, pos1, perc, m, u
            )
            .unwrap();
            writeln!(&mut bed_graph, "{}\t{}\t{}\t{}", chr, pos, pos1, perc).unwrap();
        }
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
