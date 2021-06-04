mod clip;
mod extract;
mod mbias;
mod cli;

use clap::Clap;
use bio::alphabets::dna;
use cli::ConvertToFastqOpts;
use cli::DeduplicateOpts;
use cli::ExtractPairedOpts;
use cli::ExtractSingleOpts;
use cli::Opts;
use cli::SubCommand;
use clip::{ClipperSingle, ClipperSingleConfig};
use extract::{
    write_bed_graph, write_bismark_cov, write_cytosine_report_chg, write_cytosine_report_chh,
    write_cytosine_report_cpg, Context, CytosineGenome,
};
use flate2::write::GzEncoder;
use flate2::Compression;
use mbias::{write_mbias_paired, write_mbias_single, CytosineRead};
use rust_htslib::{
    bam,
    bam::{record::Cigar, Read},
};
use std::fs::File;
use std::io::Write;
use std::{collections::HashSet, usize};


fn main() {
    let opts = Opts::parse();

    match opts.subcmd {
        SubCommand::Deduplicate(opts) => deduplicate(opts),
        SubCommand::ExtractSingle(opts) => extract_single(opts),
        SubCommand::ExtractPaired(opts) => extract_paired(opts),
        SubCommand::ConvertToFastq(opts) => convert_to_fastq(opts),
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
            start = record_2.pos() as u64 + 1;

            for c in record_1.cigar().iter() {
                match c {
                    Cigar::Match(len) | Cigar::Del(len) | Cigar::RefSkip(len) => end += *len as u64,
                    Cigar::Ins(_) | Cigar::SoftClip(_) => {}
                    _ => {
                        panic!("invalid CIGAR operation: {}", c)
                    }
                }
            }
        } else {
            end = record_2.pos() as u64;

            for c in record_2.cigar().iter() {
                match c {
                    Cigar::Match(len) | Cigar::Del(len) | Cigar::RefSkip(len) => end += *len as u64,
                    Cigar::Ins(_) | Cigar::SoftClip(_) => {}
                    _ => {
                        panic!("invalid CIGAR operation: {}", c)
                    }
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

fn extract_single(opts: ExtractSingleOpts) -> () {
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

    let config = ClipperSingleConfig {
        five_prime_clip: opts.five_prime_clip,
        three_prime_clip: opts.three_prime_clip,
        five_soft_clip: opts.five_soft_clip,
        three_soft_clip: opts.three_soft_clip,
    };
    let mut clipper = ClipperSingle::new(config);

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
            task.process_single(pos, &chr, &xm, &cigar).unwrap();
        }
    }

    if let Some(output) = opts.mbias {
        let mut writer = File::create(output).unwrap();
        write_mbias_single(cytosine_read.as_ref().unwrap(), &mut writer).unwrap();
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
    write_clipper_stats_single(&clipper);
    if let Some(task) = cytosine_read.as_ref() {
        report_read_stats(task);
    }
}

fn extract_paired(opts: ExtractPairedOpts) -> () {
    let mut input_bam = bam::Reader::from_path(&opts.input).unwrap();

    // Given BAM header, get values from key SN of tags SQ.
    let chrs = bam::Header::from_template(input_bam.header())
        .to_hashmap()
        .get("SQ")
        .unwrap()
        .iter()
        .map(|c| c.get("SN").unwrap().to_owned())
        .collect();

    let mut cytosine_read_1 = CytosineRead::new(opts.start_capacity);
    let mut cytosine_read_2 = CytosineRead::new(opts.start_capacity);

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

    // let config = ClipperPairedConfig {
    //     five_prime_clip_1: opts.five_prime_clip_1,
    //     three_prime_clip_1: opts.three_prime_clip_1,
    //     five_soft_clip_1: opts.five_soft_clip_1,
    //     three_soft_clip_1: opts.three_soft_clip_1,

    //     five_prime_clip_2: opts.five_prime_clip_2,
    //     three_prime_clip_2: opts.three_prime_clip_2,
    //     five_soft_clip_2: opts.five_soft_clip_2,
    //     three_soft_clip_2: opts.three_soft_clip_2,
    // };
    // let mut clipper = ClipperPaired::new(config);

    // Get statistics about mapped reads.
    let mut total: u32 = 0;
    let mut total_valid: u32 = 0;
    // let mut total_removed: u32 = 0;

    // let min_length = opts.min_length as usize;

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

        // Discard records that are unmapped or duplicated.
        if record_1.is_unmapped()
            || record_1.is_duplicate()
            || record_2.is_unmapped()
            || record_2.is_duplicate()
        {
            continue;
        }
        total_valid += 1;

        let chr_1 = input_bam.header().tid2name(record_1.tid() as u32);
        let chr_1 = String::from_utf8(chr_1.to_vec()).unwrap();
        let start_1 = (record_1.pos() + 1) as u64;
        let cigar_1 = record_1.cigar().to_vec();
        let xm_1 = record_1.aux(b"XM").unwrap().string().to_vec();
        let reverse_1 = is_reverse(&record_1).unwrap();

        let chr_2 = input_bam.header().tid2name(record_2.tid() as u32);
        let chr_2 = String::from_utf8(chr_2.to_vec()).unwrap();
        let start_2 = (record_2.pos() + 1) as u64;
        let cigar_2 = record_2.cigar().to_vec();
        let xm_2 = record_2.aux(b"XM").unwrap().string().to_vec();
        let reverse_2 = is_reverse(&record_2).unwrap();

        // match clipper.process(
        //     cigar_1, start_1, xm_1, reverse_1, cigar_2, start_2, xm_2, reverse_2,
        // ) {
        //     Some(((new_cigar_1, new_pos_1, new_xm_1), (new_cigar_2, new_pos_2, new_xm_2))) => {
        //         cigar_1 = new_cigar_1;
        //         start_1 = new_pos_1;
        //         xm_1 = new_xm_1;
        //         cigar_2 = new_cigar_2;
        //         start_2 = new_pos_2;
        //         xm_2 = new_xm_2;
        //     }
        //     None => {
        //         total_removed += 1;
        //         continue;
        //     }
        // }

        // if xm_1.len() < min_length || xm_2.len() < min_length {
        //     total_removed += 1;
        //     continue;
        // }

        cytosine_read_1.process(&xm_1, reverse_1, &cigar_1).unwrap();
        cytosine_read_2.process(&xm_2, reverse_2, &cigar_2).unwrap();

        if let Some(task) = cytosine_genome.as_mut() {
            task.process_paired(
                start_1, &chr_1, &xm_1, &cigar_1, start_2, &chr_2, &xm_2, &cigar_2, reverse_2,
            )
            .unwrap();
        }
    }

    if let Some(output) = opts.mbias {
        let mut writer = File::create(output).unwrap();
        write_mbias_paired(&cytosine_read_1, &cytosine_read_2, &mut writer).unwrap();
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

    // println!(
    //     "Removed:                {} ({} %)",
    //     total_removed,
    //     total_removed as f32 / total as f32 * 100.0
    // );
    // write_clipper_stats_paired(&clipper);
    report_read_stats(&cytosine_read_1);
}

fn write_clipper_stats_single(clipper: &ClipperSingle) {
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

// fn write_clipper_stats_paired(clipper: &ClipperPaired) {
//     let total = clipper.total() as f32;
//     println!(
//         "M:                      {} ({:.2} %)",
//         clipper.total_m(),
//         clipper.total_m() as f32 / total * 100.0
//     );
//     println!(
//         "SM:                     {} ({:.2} %)",
//         clipper.total_sm(),
//         clipper.total_sm() as f32 / total * 100.0
//     );
//     println!(
//         "MS:                     {} ({:.2} %)",
//         clipper.total_ms(),
//         clipper.total_ms() as f32 / total * 100.0
//     );
//     println!(
//         "SMS:                    {} ({:.2} %)",
//         clipper.total_sms(),
//         clipper.total_sms() as f32 / total * 100.0
//     );
// }

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

fn convert_to_fastq(opts: ConvertToFastqOpts) {
    // Open input and output files.
    let mut input = bam::Reader::from_path(&opts.input).unwrap();
    let mut output = GzEncoder::new(File::create(&opts.output).unwrap(), Compression::default());
    let header = bam::Header::from_template(input.header());
    let mut remove = bam::Writer::from_path(&opts.remove, &header, bam::Format::BAM).unwrap();

    // Count CIGAR string variants.
    let mut total_m: u32 = 0;
    let mut total_sm: u32 = 0;
    let mut total_ms: u32 = 0;
    let mut total_sms: u32 = 0;

    // Summary statistics.
    // Read is invalid when is unmapped or duplicated.
    // Read is removed when CIGAR string is invalid or less than threshold.
    let mut total: u32 = 0;
    let mut total_invalid: u32 = 0;
    let mut total_removed: u32 = 0;

    // Count read alignment orientation.
    let mut total_ot: u32 = 0;
    let mut total_ctot: u32 = 0;
    let mut total_ob: u32 = 0;
    let mut total_ctob: u32 = 0;

    let mut record = bam::Record::new();
    while let Some(result) = input.read(&mut record) {
        result.unwrap();
        total += 1;

        // Discard records that are unmapped, duplicated or CIGAR string not present.
        if record.is_unmapped() || record.is_duplicate() || record.cigar_len() == 0 {
            total_invalid += 1;
            continue;
        }

        // Get Bismark tags about read conversion (XR) and genome conversion (XG).
        //   XR=CT and XG=CT -> original top strand (OT)
        //   XR=GA and XG=CT -> complementary to original top strand (CTOT)
        //   XR=CT and XG=GA -> original bottom strand (OB)
        //   XR=GA and XG=GA -> complementary to original bottom strand (CTOB)
        // OB and CTOB reads are mapped as reverse
        let xr = record.aux(b"XR").expect("Missing XR tag").string();
        let xg = record.aux(b"XG").expect("Missing XG tag").string();
        let reverse = if xr == b"CT" && xg == b"CT" {
            total_ot += 1;
            false
        } else if xr == b"GA" && xg == b"GA" {
            total_ctob += 1;
            false
        } else if xr == b"CT" && xg == b"GA" {
            total_ob += 1;
            true
        } else if xr == b"GA" && xg == b"CT" {
            total_ctot += 1;
            true
        } else {
            panic!("Invalid tags: XR={:?} XG={:?}", xr, xg);
        };

        // Get CIGAR string and reverse if needed.
        let mut cigar = record.cigar().to_vec();
        if reverse {
            cigar.reverse();
        }

        // Get number of bases that were soft-clipped at both ends.
        let leading_softclips = cigar.first().map_or(0, |cigar| {
            if let Cigar::SoftClip(s) = cigar {
                *s as i32
            } else {
                0
            }
        });

        let trailing_softclips = cigar.last().map_or(0, |cigar| {
            if let Cigar::SoftClip(s) = cigar {
                *s as i32
            } else {
                0
            }
        });

        // Calculate range of valid reads according to soft-clipped bases.
        // If both leading and trailing are 0 then it is a perfect match (M).
        // If leading is greater than 0 and trailing is 0, then it is (SM).
        // If leading is 0 and tailing is greater than 0, then it is (MS).
        // Otherwise it is SMS.
        let mut start = 0;
        let mut end = record.seq_len() as i32 - 1;
        let five_soft_clip = opts.five_soft_clip as i32;
        let three_soft_clip = opts.three_soft_clip as i32;
        if leading_softclips == 0 && trailing_softclips == 0 {
            total_m += 1;
        } else if leading_softclips > 0 && trailing_softclips == 0 {
            total_sm += 1;
            start += leading_softclips + five_soft_clip;
        } else if leading_softclips == 0 && trailing_softclips > 0 {
            total_ms += 1;
            end -= trailing_softclips + three_soft_clip;
        } else {
            total_sms += 1;
            start += leading_softclips + five_soft_clip;
            end -= trailing_softclips + three_soft_clip;
        }

        // Get length of match alignment read.
        // Discard reads that have length less than threshold.
        // Write to removed BAM file.
        let length = end - start;

        if length <= 0 || length < opts.min_length as i32 {
            total_removed += 1;
            remove.write(&record).unwrap();
            continue;
        }

        // Get sequence read and quality.
        // Sum 33 to quality value to get Phred scores.
        let mut seq = record.seq().as_bytes().to_owned();
        let mut qual: Vec<u8> = record.qual().iter().map(|x| x + 33).collect();

        // Reverse read and quality.
        if reverse {
            seq = dna::revcomp(seq);
            qual.reverse();
        }

        // Clip sequence read and quality.
        // Write filtered read as FASTQ format.
        let start = start as usize;
        let end = end as usize;
        output.write_all(b"@").unwrap();
        output.write_all(record.qname()).unwrap();
        output.write_all(b"\n").unwrap();
        output.write_all(&seq[start..end]).unwrap();
        output.write_all(b"\n+\n").unwrap();
        output.write_all(&qual[start..end]).unwrap();
        output.write_all(b"\n").unwrap();
    }

    output.flush().unwrap();

    // Input parameters.
    println!("Input:        {}", opts.input.to_str().unwrap());
    println!("Output:       {}", opts.output.to_str().unwrap());
    println!("Removed:      {}", opts.remove.to_str().unwrap());
    println!("5' soft clip: {}", opts.five_soft_clip);
    println!("3' soft clip: {}", opts.three_soft_clip);
    println!("Min length:   {}", opts.min_length);

    // Summary statistics.
    let total_valid = total - total_invalid;
    let total_passed = total_valid - total_removed;
    println!("Total:        {}", total);
    let total = total as f32;
    println!(
        "Invalid:      {} ({:.2} %)",
        total_invalid,
        total_invalid as f32 / total * 100.0
    );
    println!(
        "Passed:       {} ({:.2} %)",
        total_passed,
        total_passed as f32 / total * 100.0
    );
    println!(
        "Removed:      {} ({:.2} %)",
        total_removed,
        total_removed as f32 / total * 100.0
    );
    // CIGAR string variants.
    println!(
        "M:            {} ({:.2} %)",
        total_m,
        total_m as f32 / total * 100.0
    );
    println!(
        "SM:           {} ({:.2} %)",
        total_sm,
        total_sm as f32 / total * 100.0
    );
    println!(
        "MS:           {} ({:.2} %)",
        total_ms,
        total_ms as f32 / total * 100.0
    );
    println!(
        "SMS:          {} ({:.2} %)",
        total_sms,
        total_sms as f32 / total * 100.0
    );

    // Read alignment orientation.
    let total_valid = total_valid as f32;
    println!(
        "OT:           {} ({:.2} %)",
        total_ot,
        total_ot as f32 / total_valid * 100.0
    );
    println!(
        "CTOB:         {} ({:.2} %)",
        total_ctob,
        total_ctob as f32 / total_valid * 100.0
    );
    println!(
        "OB:           {} ({:.2} %)",
        total_ob,
        total_ob as f32 / total_valid * 100.0
    );
    println!(
        "CTOT:         {} ({:.2} %)",
        total_ctot,
        total_ctot as f32 / total_valid * 100.0
    );
}
