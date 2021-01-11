use clap::Clap;
use flate2::write::GzEncoder;
use flate2::Compression;
use methylcx::{CytosineGenome, CytosineRead};
use rust_htslib::{bam, bam::Read};
use std::io::Write;
use std::path::PathBuf;
use std::usize;
use std::{fs::File, io};

#[derive(Clap)]
struct Opts {
    /// Input mapped reads by Bismark aligner (BAM/SAM/CRAM).
    #[clap(parse(from_os_str))]
    input: PathBuf,
    /// Counts of (un)methylated and coverage across cycles (CSV).
    #[clap(long, parse(from_os_str))]
    mbias: Option<PathBuf>,
    /// Counts of (un)methylated and coverage across genome (cov, gzipped).
    #[clap(long, parse(from_os_str))]
    bismark_cov: Option<PathBuf>,
    /// Counts of (un)methylated and coverage across genome (bedGraph, gzipped).
    #[clap(long, parse(from_os_str))]
    bed_graph: Option<PathBuf>,
    /// Starting vector capacity to store per sequencing cycles data.
    #[clap(long, default_value = "200")]
    start_capacity: usize,
}

fn main() {
    let opts: Opts = Opts::parse();

    let mut bam = bam::Reader::from_path(&opts.input).unwrap();

    // Given BAM header, get values from key SN of tags SQ.
    let chrs = bam::Header::from_template(bam.header())
        .to_hashmap()
        .get("SQ")
        .unwrap()
        .iter()
        .map(|c| c.get("SN").unwrap().to_owned())
        .collect();
    let mut cytosine_genome = CytosineGenome::new(chrs);
    let mut cytosine_read = CytosineRead::new(opts.start_capacity);

    // Get statistics about mapped reads.
    let mut total = 0;
    let mut total_valid = 0;

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

        let chr = bam.header().tid2name(record.tid() as u32).to_owned();
        let chr = String::from_utf8(chr).unwrap();
        let pos = record.pos() as u64;
        let cigar = record.cigar().to_vec();
        let xm = record.aux(b"XM").unwrap().string();
        let reverse = is_reverse(&record).unwrap();

        cytosine_read.process(xm, reverse, &cigar);
        cytosine_genome.process(pos, chr, xm, &cigar).unwrap();
    }

    if let Some(output) = opts.mbias {
        let mut writer = File::create(output).unwrap();
        write_mbias(&cytosine_read, &mut writer).unwrap();
    }

    if let Some(output) = opts.bismark_cov {
        let mut writer = GzEncoder::new(File::create(output).unwrap(), Compression::default());
        write_bed_graph(&cytosine_genome, &mut writer).unwrap();
    }

    if let Some(output) = opts.bed_graph {
        let mut writer = GzEncoder::new(File::create(output).unwrap(), Compression::default());
        write_bismark_cov(&cytosine_genome, &mut writer).unwrap();
    }

    let percent_valid = total_valid as f32 / total as f32 * 100.0;
    println!("Total:                  {}", total);
    println!(
        "Valid:                  {} ({} %)",
        total_valid, percent_valid
    );
    report_read_stats(&cytosine_read);
}

// Print result tables in CSV format to stdout.
// Tables are merged as long format and discriminated by `Context` column.
fn write_mbias(cytosine_read: &CytosineRead, writer: &mut File) -> io::Result<()> {
    writer.write_all(b"Context,Cycle,Methylated,Unmethylated,Coverage\n")?;
    for (i, (m, um)) in cytosine_read
        .cpg_m()
        .iter()
        .zip(cytosine_read.cpg_u())
        .enumerate()
    {
        writeln!(writer, "CpG,{},{},{},{}", i + 1, m, um, m + um)?;
    }
    for (i, (m, um)) in cytosine_read
        .chg_m()
        .iter()
        .zip(cytosine_read.chg_u())
        .enumerate()
    {
        writeln!(writer, "CHG,{},{},{},{}", i + 1, m, um, m + um)?;
    }
    for (i, (m, um)) in cytosine_read
        .chh_m()
        .iter()
        .zip(cytosine_read.chh_u())
        .enumerate()
    {
        writeln!(writer, "CHH,{},{},{},{}", i + 1, m, um, m + um)?;
    }
    Ok(())
}

fn write_bed_graph(
    cytosine_genome: &CytosineGenome,
    writer: &mut GzEncoder<File>,
) -> io::Result<()> {
    writer.write_all(b"track type=bedGraph\n")?;
    for (chr, xs) in cytosine_genome.cpg() {
        for (pos, (m, _, cov)) in xs {
            let perc = *m as f64 / *cov as f64 * 100.0;
            let pos1 = pos + 1;
            writeln!(writer, "{}\t{}\t{}\t{}", chr, pos, pos1, perc)?;
        }
    }
    Ok(())
}

fn write_bismark_cov(
    cytosine_genome: &CytosineGenome,
    writer: &mut GzEncoder<File>,
) -> io::Result<()> {
    for (chr, xs) in cytosine_genome.cpg() {
        for (pos, (m, u, cov)) in xs {
            let perc = *m as f64 / *cov as f64 * 100.0;
            let pos1 = pos + 1;
            writeln!(
                writer,
                "{}\t{}\t{}\t{}\t{}\t{}",
                chr, pos1, pos1, perc, m, u
            )?;
        }
    }
    return Ok(());
}

// Get Bismark tags about read conversion (XR) and genome conversion (XG)
// XR=CT and XG=CT -> original top strand (OT)
// XR=GA and XG=CT -> complementary to original top strand (CTOT)
// XR=CT and XG=GA -> original bottom strand (OB)
// XR=GA and XG=GA -> complementary to original bottom strand (CTOB)
// TODO: explain why reverse read when CT/GA or GA/CT
fn is_reverse(record: &bam::Record) -> Result<bool, String> {
    let xr = match record.aux(b"XR") {
        Some(value) => value.string(),
        None => return Err(String::from("missing XR tag")),
    };
    let xg = match record.aux(b"XG") {
        Some(value) => value.string(),
        None => return Err(String::from("missing XG tag")),
    };

    Ok((xr == b"CT" && xg == b"GA") || (xr == b"GA" && xg == b"CT"))
}

fn report_read_stats(cr: &CytosineRead) {
    let total_c = cr.total_c();
    let total_cpg = cr.total_cpg();
    let total_chg = cr.total_chg();
    let total_chh = cr.total_chh();
    let total_cpg_m = cr.total_cpg_methylated();
    let total_cpg_u = cr.total_cpg_unmethylated();
    let total_chg_m = cr.total_chg_methylated();
    let total_chg_u = cr.total_chg_unmethylated();
    let total_chh_m = cr.total_chh_methylated();
    let total_chh_u = cr.total_chh_unmethylated();
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
        total_cpg_u,
        total_cpg_u as f32 / total_cpg as f32 * 100.0
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
        total_chg_u,
        total_chg_u as f32 / total_chg as f32 * 100.0
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
        total_chh_u,
        total_chh_u as f32 / total_chh as f32 as f32 * 100.0
    );
}
