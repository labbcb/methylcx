use bio::{alphabets::dna, io::fasta};
use clap::Clap;
use flate2::read::GzDecoder;
use flate2::write::GzEncoder;
use flate2::Compression;
use methylcx::{Clipper, ClipperConfig, CytosineGenome, CytosineRead};
use rust_htslib::{
    bam,
    bam::{record::Cigar, Read},
};
use std::io::Write;
use std::path::PathBuf;
use std::{collections::HashSet, usize};
use std::{fs::File, io};

#[derive(Clap)]
struct Opts {
    #[clap(subcommand)]
    subcmd: SubCommand,
}

#[derive(Clap)]
enum SubCommand {
    Deduplicate(DeduplicateOpts),
    Extract(ExtractOpts),
}

/// Removes duplicated reads while keeping the first one
#[derive(Clap)]
struct DeduplicateOpts {
    /// Input mapped reads by Bismark aligner (BAM/SAM).
    #[clap(parse(from_os_str))]
    input: PathBuf,

    /// Output deduplicated reads (BAM).
    #[clap(parse(from_os_str))]
    output: PathBuf,

    /// Either is paired-end sequencing or single-end (default)
    #[clap(long)]
    paired: bool,
}

/// Calculates DNA methylation across genome
#[derive(Clap)]
struct ExtractOpts {
    /// Input mapped reads by Bismark aligner (BAM/SAM).
    #[clap(parse(from_os_str))]
    input: PathBuf,
    /// Clip first bases of aligned sequence in end-to-end mode (5' orientation).
    #[clap(long, default_value = "0")]
    five_prime_clip: u32,
    /// Clip last bases of aligned sequence in end-to-end mode (3' orientation).
    #[clap(long, default_value = "0")]
    three_prime_clip: u32,
    /// Clip first bases of aligned sequence in local mode (5' orientation).
    #[clap(long, default_value = "0")]
    five_soft_clip: u32,
    /// Clip last bases of aligned sequence in local mode (3' orientation).
    #[clap(long, default_value = "0")]
    three_soft_clip: u32,
    /// Minimum read length.
    #[clap(long, default_value = "0")]
    min_length: u32,
    /// Counts of (un)methylated and coverage across cycles (CSV).
    #[clap(long, parse(from_os_str))]
    mbias: Option<PathBuf>,

    /// Counts of (un)methylated and coverage across genome in CpG-context (cov, gzipped).
    #[clap(long, parse(from_os_str))]
    cpg_bismark_cov: Option<PathBuf>,
    /// Counts of (un)methylated and coverage across genome in CHG-context (cov, gzipped).
    #[clap(long, parse(from_os_str))]
    chg_bismark_cov: Option<PathBuf>,
    /// Counts of (un)methylated and coverage across genome in CHH-context (cov, gzipped).
    #[clap(long, parse(from_os_str))]
    chh_bismark_cov: Option<PathBuf>,

    /// Counts of (un)methylated and coverage across genome in CpG-context (bedGraph, gzipped).
    #[clap(long, parse(from_os_str))]
    cpg_bed_graph: Option<PathBuf>,
    /// Counts of (un)methylated and coverage across genome in CHG-context (bedGraph, gzipped).
    #[clap(long, parse(from_os_str))]
    chg_bed_graph: Option<PathBuf>,
    /// Counts of (un)methylated and coverage across genome in CHH-context (bedGraph, gzipped).
    #[clap(long, parse(from_os_str))]
    chh_bed_graph: Option<PathBuf>,

    /// Minimum coverage to report loci. To output every loci see --cytosine-report.
    #[clap(long, default_value = "1")]
    min_coverage: u32,

    /// Strand-specific counts of (un)methylated CpG across genome (cytosine-report, gzipped).
    #[clap(long, parse(from_os_str))]
    cpg_cytosine_report: Option<PathBuf>,
    /// Strand-specific counts of (un)methylated CHG across genome (cytosine-report, gzipped).
    #[clap(long, parse(from_os_str))]
    chg_cytosine_report: Option<PathBuf>,
    /// Strand-specific counts of (un)methylated CHH across genome (cytosine-report, gzipped).
    #[clap(long, parse(from_os_str))]
    chh_cytosine_report: Option<PathBuf>,

    /// Reference genome file (FASTA format, gzipped). Only required for --cytosine-report.
    #[clap(long, parse(from_os_str))]
    genome: Option<PathBuf>,

    /// Starting vector capacity to store per sequencing cycles data.
    #[clap(long, default_value = "200")]
    start_capacity: usize,
}

fn main() {
    let opts: Opts = Opts::parse();

    match opts.subcmd {
        SubCommand::Deduplicate(opts) => deduplicate(opts),
        SubCommand::Extract(opts) => extract(opts),
    }
}

fn deduplicate(opts: DeduplicateOpts) {
    let mut input_bam = bam::Reader::from_path(&opts.input).unwrap();
    let header = bam::Header::from_template(input_bam.header());

    let mut output_bam = bam::Writer::from_path(opts.output, &header, bam::Format::BAM).unwrap();

    let mut unique_seqs = HashSet::new();

    let total = if opts.paired {
        deduplicate_paired(&mut input_bam, &mut unique_seqs, &mut output_bam)
    } else {
        deduplicate_single(&mut input_bam, &mut unique_seqs, &mut output_bam)
    };

    let total_unique = unique_seqs.len() as u32;
    let total_duplicated = total - total_unique;

    eprintln!("Total: {}", total);
    let total = total as f32;
    eprintln!(
        "Unique: {} ({:.2} %)",
        total_unique,
        total_unique as f32 / total * 100.0
    );
    eprintln!(
        "Duplicates: {} ({:.2} %)",
        total_duplicated,
        total_duplicated as f32 / total * 100.0
    );
}

fn deduplicate_paired(
    input_bam: &mut bam::Reader,
    unique_seqs: &mut HashSet<String>,
    output_bam: &mut bam::Writer,
) -> u32 {
    let mut total = 0;

    let mut record_1 = bam::Record::new();
    let mut record_2 = bam::Record::new();
    loop {
        let result_1 = input_bam.read(&mut record_1);
        if let Some(result) = result_1 {
            result.unwrap();
        } else {
            break;
        }

        let result_2 = input_bam.read(&mut record_2);
        if let Some(result) = result_2 {
            result.unwrap();
        } else {
            panic!(
                "missing read mate for {}",
                std::str::from_utf8(record_1.qname()).unwrap()
            );
        }

        total += 1;

        let chr = input_bam.header().tid2name(record_1.tid() as u32);
        let chr = String::from_utf8(chr.to_vec()).unwrap();
        let mut start = record_1.pos() as u64 + 1;
        let strand = get_strand(&record_1).unwrap();

        let mut end: u64;
        if strand == "CTOT" || strand == "OB" {
            end = start - 1;
            start = record_2.pos() as u64;
        } else {
            end = record_2.pos() as u64 - 1;
        }

        for c in record_2.cigar().iter() {
            match c {
                Cigar::Match(len) | Cigar::Del(len) | Cigar::RefSkip(len) => end += *len as u64,
                Cigar::Ins(_) | Cigar::SoftClip(_) => {}
                _ => {
                    panic!("invalid CIGAR operation: {}", c)
                }
            }
        }

        let id = format!("{}:{}:{}:{}", strand, chr, start, end);
        if !unique_seqs.contains(&id) {
            unique_seqs.insert(id);

            output_bam.write(&record_1).unwrap();
            output_bam.write(&record_2).unwrap();
        }
    }

    total
}

fn deduplicate_single(
    input_bam: &mut bam::Reader,
    unique_seqs: &mut HashSet<String>,
    output_bam: &mut bam::Writer,
) -> u32 {
    let mut total = 0;

    let mut record = bam::Record::new();
    while let Some(result) = input_bam.read(&mut record) {
        result.unwrap();

        total += 1;

        let chr = input_bam.header().tid2name(record.tid() as u32);
        let chr = std::str::from_utf8(chr).unwrap();
        let mut pos = record.pos() as u32 + 1;
        let strand = get_strand(&record).unwrap();

        if strand == "CTOT" || strand == "OB" {
            pos -= 1;
            for c in record.cigar().iter() {
                match c {
                    Cigar::Match(len) | Cigar::Del(len) | Cigar::RefSkip(len) => pos += *len,
                    Cigar::Ins(_) | Cigar::SoftClip(_) => {}
                    _ => {
                        panic!("invalid CIGAR operation: {}", c)
                    }
                }
            }
        }

        let id = format!("{}:{}:{}", strand, chr, pos);
        if !unique_seqs.contains(&id) {
            unique_seqs.insert(id);

            output_bam.write(&record).unwrap();
        }
    }

    total
}

fn extract(opts: ExtractOpts) -> () {
    let mut bam = bam::Reader::from_path(&opts.input).unwrap();

    // Given BAM header, get values from key SN of tags SQ.
    let chrs = bam::Header::from_template(bam.header())
        .to_hashmap()
        .get("SQ")
        .unwrap()
        .iter()
        .map(|c| c.get("SN").unwrap().to_owned())
        .collect();

    // Instanciate processing tasks on demand.
    let mut cytosine_read = if opts.mbias.is_some() {
        Some(CytosineRead::new(opts.start_capacity))
    } else {
        None
    };

    let run_cytosine_genome = opts.cpg_bismark_cov.is_some()
        || opts.chg_bismark_cov.is_some()
        || opts.chh_bismark_cov.is_some()
        || opts.cpg_bed_graph.is_some()
        || opts.chg_bed_graph.is_some()
        || opts.chh_bed_graph.is_some()
        || opts.cpg_cytosine_report.is_some()
        || opts.chg_cytosine_report.is_some()
        || opts.chh_cytosine_report.is_some();
    let mut cytosine_genome = if run_cytosine_genome {
        if opts.cpg_cytosine_report.is_some() && opts.genome.is_none() {
            panic!("mising reference genome file");
        }
        Some(CytosineGenome::new(chrs))
    } else {
        None
    };

    let config = ClipperConfig {
        five_prime_clip: opts.five_prime_clip,
        three_prime_clip: opts.three_prime_clip,
        five_soft_clip: opts.five_soft_clip,
        three_soft_clip: opts.three_soft_clip,
    };
    let mut clipper = Clipper::new(config);

    // Get statistics about mapped reads.
    let mut total: u32 = 0;
    let mut total_valid: u32 = 0;
    let mut total_removed: u32 = 0;

    let min_length = opts.min_length as usize;

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

        let chr = bam.header().tid2name(record.tid() as u32);
        let chr = String::from_utf8(chr.to_vec()).unwrap();

        let mut pos = (record.pos() + 1) as u64;
        let mut cigar = record.cigar().to_vec();
        let mut xm = record.aux(b"XM").unwrap().string().to_vec();
        let reverse = is_reverse(&record).unwrap();

        match clipper.process(cigar, pos, xm, reverse) {
            Some((new_cigar, new_pos, new_xm)) => {
                cigar = new_cigar;
                pos = new_pos;
                xm = new_xm;
            }
            None => {
                total_removed += 1;
                continue;
            }
        }

        if xm.len() < min_length {
            total_removed += 1;
            continue;
        }

        if let Some(task) = cytosine_read.as_mut() {
            task.process(&xm, reverse, &cigar).unwrap();
        }
        if let Some(task) = cytosine_genome.as_mut() {
            task.process(pos, &chr, &xm, &cigar).unwrap();
        }
    }

    if let Some(output) = opts.mbias {
        let mut writer = File::create(output).unwrap();
        write_mbias(cytosine_read.as_ref().unwrap(), &mut writer).unwrap();
    }

    if let Some(output) = opts.cpg_bismark_cov {
        let mut writer = GzEncoder::new(File::create(output).unwrap(), Compression::default());
        write_bismark_cov(
            cytosine_genome.as_ref().unwrap(),
            opts.min_coverage,
            Context::CpG,
            &mut writer,
        )
        .unwrap();
    }

    if let Some(output) = opts.chg_bismark_cov {
        let mut writer = GzEncoder::new(File::create(output).unwrap(), Compression::default());
        write_bismark_cov(
            cytosine_genome.as_ref().unwrap(),
            opts.min_coverage,
            Context::CHG,
            &mut writer,
        )
        .unwrap();
    }

    if let Some(output) = opts.chh_bismark_cov {
        let mut writer = GzEncoder::new(File::create(output).unwrap(), Compression::default());
        write_bismark_cov(
            cytosine_genome.as_ref().unwrap(),
            opts.min_coverage,
            Context::CHH,
            &mut writer,
        )
        .unwrap();
    }

    if let Some(output) = opts.cpg_bed_graph {
        let mut writer = GzEncoder::new(File::create(output).unwrap(), Compression::default());
        write_bed_graph(
            cytosine_genome.as_ref().unwrap(),
            opts.min_coverage,
            Context::CpG,
            &mut writer,
        )
        .unwrap();
    }

    if let Some(output) = opts.chg_bed_graph {
        let mut writer = GzEncoder::new(File::create(output).unwrap(), Compression::default());
        write_bed_graph(
            cytosine_genome.as_ref().unwrap(),
            opts.min_coverage,
            Context::CHG,
            &mut writer,
        )
        .unwrap();
    }

    if let Some(output) = opts.chh_bed_graph {
        let mut writer = GzEncoder::new(File::create(output).unwrap(), Compression::default());
        write_bed_graph(
            cytosine_genome.as_ref().unwrap(),
            opts.min_coverage,
            Context::CHH,
            &mut writer,
        )
        .unwrap();
    }

    if let Some(output) = opts.cpg_cytosine_report {
        let genome = File::open(opts.genome.as_ref().unwrap()).unwrap();
        let mut writer = GzEncoder::new(File::create(output).unwrap(), Compression::default());
        write_cytosine_report_cpg(cytosine_genome.as_ref().unwrap(), &genome, &mut writer).unwrap();
    }

    if let Some(output) = opts.chg_cytosine_report {
        let genome = File::open(opts.genome.as_ref().unwrap()).unwrap();
        let mut writer = GzEncoder::new(File::create(output).unwrap(), Compression::default());
        write_cytosine_report_chg(cytosine_genome.as_ref().unwrap(), &genome, &mut writer).unwrap();
    }

    if let Some(output) = opts.chh_cytosine_report {
        let genome = File::open(opts.genome.as_ref().unwrap()).unwrap();
        let mut writer = GzEncoder::new(File::create(output).unwrap(), Compression::default());
        write_cytosine_report_chh(cytosine_genome.as_ref().unwrap(), &genome, &mut writer).unwrap();
    }

    let percent_valid = total_valid as f32 / total as f32 * 100.0;
    println!("Total:                  {}", total);
    println!(
        "Valid:                  {} ({} %)",
        total_valid, percent_valid
    );

    println!(
        "Removed:                {} ({} %)",
        total_removed,
        total_removed as f32 / total as f32 * 100.0
    );
    write_clipper_stats(&clipper);
    if let Some(task) = cytosine_read.as_ref() {
        report_read_stats(task);
    }
}

#[derive(PartialEq)]
enum Context {
    CpG,
    CHG,
    CHH,
}

impl Context {
    fn parse(triplet: &[u8]) -> Option<Context> {
        if triplet[0] != b'C' && triplet[0] != b'c' {
            None
        } else if triplet[1] == b'G' || triplet[1] == b'g' {
            Some(Context::CpG)
        } else if triplet[2] == b'G' || triplet[2] == b'g' {
            Some(Context::CHG)
        } else {
            Some(Context::CHH)
        }
    }

    // TODO: could be `Context::parse(triplet).contains(context)`?
    fn is(triplet: &[u8], context: Context) -> bool {
        match Context::parse(triplet) {
            Some(ctx) => ctx == context,
            None => false,
        }
    }
}

fn write_clipper_stats(clipper: &Clipper) {
    let total = clipper.total() as f32;
    println!(
        "M:                      {} ({:.2} %)",
        clipper.total_m(),
        clipper.total_m() as f32 / total * 100.0
    );
    println!(
        "SM:                     {} ({:.2} %)",
        clipper.total_sm(),
        clipper.total_sm() as f32 / total * 100.0
    );
    println!(
        "MS:                     {} ({:.2} %)",
        clipper.total_ms(),
        clipper.total_ms() as f32 / total * 100.0
    );
    println!(
        "SMS:                    {} ({:.2} %)",
        clipper.total_sms(),
        clipper.total_sms() as f32 / total * 100.0
    );
}

fn write_cytosine_report_cpg(
    cytosine_genome: &CytosineGenome,
    genome: &File,
    writer: &mut GzEncoder<File>,
) -> io::Result<()> {
    let reader = fasta::Reader::new(GzDecoder::new(genome));

    let mut records = reader.records();
    while let Some(Ok(record)) = records.next() {
        let id = record.id();
        let seq = record.seq();

        let cpg = match cytosine_genome.cpg().get(id) {
            Some(value) => value,
            None => {
                continue;
            }
        };

        let mut pos = 0;
        for idx in 0..seq.len() - 2 {
            pos += 1;

            if !((seq[idx] == b'C' || seq[idx] == b'c')
                && (seq[idx + 1] == b'G' || seq[idx + 1] == b'g'))
            {
                continue;
            }

            let (m, u, _) = cpg.get(&pos).unwrap_or(&(0, 0, 0));

            writeln!(
                writer,
                "{}\t{}\t+\t{}\t{}\tCG\t{}",
                id,
                pos,
                m,
                u,
                std::str::from_utf8(&seq[idx..idx + 3]).unwrap()
            )?;

            let (m, u, _) = cpg.get(&(pos + 1)).unwrap_or(&(0, 0, 0));

            writeln!(
                writer,
                "{}\t{}\t-\t{}\t{}\tCG\t{}",
                id,
                pos + 1,
                m,
                u,
                std::str::from_utf8(&dna::revcomp(&seq[idx - 1..idx + 2])).unwrap()
            )?;
        }
    }

    Ok(())
}

fn write_cytosine_report_chg(
    cytosine_genome: &CytosineGenome,
    genome: &File,
    writer: &mut GzEncoder<File>,
) -> io::Result<()> {
    let reader = fasta::Reader::new(GzDecoder::new(genome));
    let mut records = reader.records();
    while let Some(Ok(record)) = records.next() {
        let id = record.id();
        let seq = record.seq();

        let chg = match cytosine_genome.chg().get(id) {
            Some(value) => value,
            None => {
                continue;
            }
        };

        let mut pos = 0;
        for idx in 0..seq.len() - 2 {
            pos += 1;

            let triplet = &seq[idx..idx + 3];
            if !Context::is(triplet, Context::CHG) {
                continue;
            }

            let (m, u, _) = chg.get(&pos).unwrap_or(&(0, 0, 0));

            writeln!(
                writer,
                "{}\t{}\t+\t{}\t{}\tCG\t{}",
                id,
                pos,
                m,
                u,
                std::str::from_utf8(triplet).unwrap()
            )?;

            let (m, u, _) = chg.get(&(pos + 2)).unwrap_or(&(0, 0, 0));

            writeln!(
                writer,
                "{}\t{}\t-\t{}\t{}\tCG\t{}",
                id,
                pos + 1,
                m,
                u,
                std::str::from_utf8(&dna::revcomp(triplet)).unwrap()
            )?;
        }
    }

    Ok(())
}

fn write_cytosine_report_chh(
    cytosine_genome: &CytosineGenome,
    genome: &File,
    writer: &mut GzEncoder<File>,
) -> io::Result<()> {
    let reader = fasta::Reader::new(GzDecoder::new(genome));
    let mut records = reader.records();
    while let Some(Ok(record)) = records.next() {
        let id = record.id();
        let seq = record.seq();

        let chh = match cytosine_genome.chh().get(id) {
            Some(value) => value,
            None => {
                continue;
            }
        };

        let mut pos = 0;
        for triplet in seq.windows(3) {
            pos += 1;

            let triplet_rev = dna::revcomp(triplet);
            let (pos, strand, triplet) = if Context::is(triplet, Context::CHH) {
                (pos, "+", triplet)
            } else if Context::is(&triplet_rev, Context::CHH) {
                (pos + 2, "-", triplet_rev.as_slice())
            } else {
                continue;
            };

            let (m, u, _) = chh.get(&pos).unwrap_or(&(0, 0, 0));

            writeln!(
                writer,
                "{}\t{}\t{}\t{}\t{}\tCHH\t{}",
                id,
                pos,
                strand,
                m,
                u,
                std::str::from_utf8(triplet).unwrap()
            )?;
        }
    }

    Ok(())
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
    min_coverage: u32,
    context: Context,
    writer: &mut GzEncoder<File>,
) -> io::Result<()> {
    assert_ne!(min_coverage, 0, "min_coverage must be larger than 0");

    let map = match context {
        Context::CpG => cytosine_genome.cpg(),
        Context::CHG => cytosine_genome.chg(),
        Context::CHH => cytosine_genome.chh(),
    };
    writer.write_all(b"track type=bedGraph\n")?;
    for chr in cytosine_genome.chrs() {
        let xs = map.get(chr).unwrap();
        for (pos, (m, _, cov)) in xs {
            if *cov < min_coverage {
                continue;
            }
            let perc = *m as f64 / *cov as f64 * 100.0;
            writeln!(writer, "{}\t{}\t{}\t{}", chr, pos - 1, pos, perc)?;
        }
    }
    Ok(())
}

fn write_bismark_cov(
    cytosine_genome: &CytosineGenome,
    min_coverage: u32,
    context: Context,
    writer: &mut GzEncoder<File>,
) -> io::Result<()> {
    assert_ne!(min_coverage, 0, "min_coverage must be larger than 0");

    let map = match context {
        Context::CpG => cytosine_genome.cpg(),
        Context::CHG => cytosine_genome.chg(),
        Context::CHH => cytosine_genome.chh(),
    };
    for chr in cytosine_genome.chrs() {
        let xs = map.get(chr).unwrap();
        for (pos, (m, u, cov)) in xs {
            if *cov < min_coverage {
                continue;
            }
            let perc = *m as f64 / *cov as f64 * 100.0;
            writeln!(writer, "{}\t{}\t{}\t{}\t{}\t{}", chr, pos, pos, perc, m, u)?;
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

fn get_strand(record: &bam::Record) -> Result<&'static str, String> {
    let xr = match record.aux(b"XR") {
        Some(value) => value.string(),
        None => return Err(String::from("missing XR tag")),
    };
    let xg = match record.aux(b"XG") {
        Some(value) => value.string(),
        None => return Err(String::from("missing XG tag")),
    };

    return if xr == b"CT" && xg == b"CT" {
        Ok("OT")
    } else if xr == b"GA" && xg == b"CT" {
        Ok("CTOT")
    } else if xr == b"GA" && xg == b"GA" {
        Ok("CTOB")
    } else if xr == b"CT" && xg == b"GA" {
        Ok("OB")
    } else {
        Err(format!(
            "invalid XR={} or XG={}",
            std::str::from_utf8(xr).unwrap(),
            std::str::from_utf8(xg).unwrap()
        ))
    };
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
