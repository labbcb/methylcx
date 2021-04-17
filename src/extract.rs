use std::{collections::BTreeMap, usize};

use bio::{alphabets::dna, io::fasta};
use flate2::read::GzDecoder;
use flate2::write::GzEncoder;
use rust_htslib::bam::record::Cigar;
use std::io::Write;
use std::{fs::File, io};

#[derive(PartialEq)]
pub enum Context {
    CpG,
    CHG,
    CHH,
}

impl Context {
    pub fn parse(triplet: &[u8]) -> Option<Context> {
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
    pub fn is(triplet: &[u8], context: Context) -> bool {
        match Context::parse(triplet) {
            Some(ctx) => ctx == context,
            None => false,
        }
    }
}
pub struct CytosineGenome {
    chrs: Vec<String>,
    cpg: BTreeMap<String, BTreeMap<u64, (u32, u32, u32)>>,
    chg: BTreeMap<String, BTreeMap<u64, (u32, u32, u32)>>,
    chh: BTreeMap<String, BTreeMap<u64, (u32, u32, u32)>>,
    unknown: BTreeMap<String, BTreeMap<u64, (u32, u32, u32)>>,
}

impl CytosineGenome {
    pub fn new(chrs: Vec<String>) -> Self {
        let mut cpg = BTreeMap::new();
        let mut chg = BTreeMap::new();
        let mut chh = BTreeMap::new();
        let mut unknown = BTreeMap::new();

        for chr in &chrs {
            cpg.insert(chr.to_string(), BTreeMap::new());
            chg.insert(chr.to_string(), BTreeMap::new());
            chh.insert(chr.to_string(), BTreeMap::new());
            unknown.insert(chr.to_string(), BTreeMap::new());
        }

        CytosineGenome {
            chrs,
            cpg,
            chg,
            chh,
            unknown,
        }
    }

    pub fn process_single(
        &mut self,
        pos: u64,
        chr: &String,
        xm: &[u8],
        cigar: &[Cigar],
    ) -> Result<(), String> {
        if !self.chrs.contains(chr) {
            return Err(format!("invalid chromosome {}", chr));
        }

        let cpg = self.cpg.get_mut(chr).unwrap();
        let chg = self.chg.get_mut(chr).unwrap();
        let chh = self.chh.get_mut(chr).unwrap();
        let unknown = self.unknown.get_mut(chr).unwrap();

        let mut pos = pos;
        let mut xm_iter = xm.iter();
        for c in cigar.iter() {
            match c {
                Cigar::Match(len) => {
                    let len = *len as usize;
                    for _ in 0..len {
                        let &b = xm_iter.next().unwrap();
                        if b == b'.' {
                            pos += 1;
                            continue;
                        }

                        let (m, u, c) = match b {
                            b'X' | b'x' => chg.entry(pos).or_insert((0, 0, 0)),
                            b'H' | b'h' => chh.entry(pos).or_insert((0, 0, 0)),
                            b'Z' | b'z' => cpg.entry(pos).or_insert((0, 0, 0)),
                            b'U' | b'u' => unknown.entry(pos).or_insert((0, 0, 0)),
                            _ => return Err(format!("invalid XM key {}", b as char)),
                        };

                        if matches!(b, b'X' | b'H' | b'Z' | b'U') {
                            *m += 1;
                        } else {
                            *u += 1;
                        }
                        *c += 1;
                        pos += 1;
                    }
                }
                Cigar::Ins(len) | Cigar::SoftClip(len) => {
                    for _ in 0..*len as usize {
                        xm_iter.next();
                    }
                }
                Cigar::Del(len) => {
                    pos += *len as u64;
                }
                _ => return Err(format!("invalid CIGAR operation {}", c.char())),
            }
        }
        Ok(())
    }

    pub fn process_paired(
        &mut self,
        start_1: u64,
        chr_1: &String,
        xm_1: &[u8],
        cigar_1: &[Cigar],
        start_2: u64,
        chr_2: &String,
        xm_2: &[u8],
        cigar_2: &[Cigar],
        reverse_2: bool,
    ) -> Result<(), String> {
        if !self.chrs.contains(chr_1) {
            return Err(format!("invalid chromosome '{}'", chr_1));
        }

        let cpg = self.cpg.get_mut(chr_1).unwrap();
        let chg = self.chg.get_mut(chr_1).unwrap();
        let chh = self.chh.get_mut(chr_1).unwrap();
        let unknown = self.unknown.get_mut(chr_1).unwrap();

        let mut pos = start_1;
        let mut xm_iter = xm_1.iter();
        for c in cigar_1.iter() {
            match c {
                Cigar::Match(len) => {
                    let len = *len as usize;
                    for _ in 0..len {
                        let &b = xm_iter.next().unwrap();
                        if b == b'.' {
                            pos += 1;
                            continue;
                        }

                        let (m, u, c) = match b {
                            b'X' | b'x' => chg.entry(pos).or_insert((0, 0, 0)),
                            b'H' | b'h' => chh.entry(pos).or_insert((0, 0, 0)),
                            b'Z' | b'z' => cpg.entry(pos).or_insert((0, 0, 0)),
                            b'U' | b'u' => unknown.entry(pos).or_insert((0, 0, 0)),
                            _ => return Err(format!("invalid XM key {}", b as char)),
                        };

                        if matches!(b, b'X' | b'H' | b'Z' | b'U') {
                            *m += 1;
                        } else {
                            *u += 1;
                        }
                        *c += 1;
                        pos += 1;
                    }
                }
                Cigar::Ins(len) | Cigar::SoftClip(len) => {
                    for _ in 0..*len as usize {
                        xm_iter.next();
                    }
                }
                Cigar::Del(len) => {
                    pos += *len as u64;
                }
                _ => return Err(format!("invalid CIGAR operation {}", c.char())),
            }
        }

        let cpg = self.cpg.get_mut(chr_2).unwrap();
        let chg = self.chg.get_mut(chr_2).unwrap();
        let chh = self.chh.get_mut(chr_2).unwrap();
        let unknown = self.unknown.get_mut(chr_2).unwrap();

        if reverse_2 {
            let end_1 = pos - 1;
            let mut pos = end_cigar(start_2, cigar_2);
            let mut xm_iter = xm_2.iter().rev();
            for c in cigar_2.iter().rev() {
                match c {
                    Cigar::Match(len) => {
                        let len = *len as usize;
                        for _ in 0..len {
                            let &b = xm_iter.next().unwrap();
                            if pos <= end_1 {
                                break;
                            }
                            if b == b'.' {
                                pos -= 1;
                                continue;
                            }

                            let (m, u, c) = match b {
                                b'X' | b'x' => chg.entry(pos).or_insert((0, 0, 0)),
                                b'H' | b'h' => chh.entry(pos).or_insert((0, 0, 0)),
                                b'Z' | b'z' => cpg.entry(pos).or_insert((0, 0, 0)),
                                b'U' | b'u' => unknown.entry(pos).or_insert((0, 0, 0)),
                                _ => return Err(format!("invalid XM key {}", b as char)),
                            };

                            if matches!(b, b'X' | b'H' | b'Z' | b'U') {
                                *m += 1;
                            } else {
                                *u += 1;
                            }
                            *c += 1;
                            pos -= 1;
                        }
                    }
                    Cigar::Ins(len) | Cigar::SoftClip(len) => {
                        for _ in 0..*len as usize {
                            xm_iter.next();
                        }
                    }
                    Cigar::Del(len) => {
                        pos -= *len as u64;
                    }
                    _ => return Err(format!("invalid CIGAR operation {}", c.char())),
                }
            }
        } else {
            let mut pos = start_2;
            let mut xm_iter = xm_2.iter();
            for c in cigar_2.iter() {
                match c {
                    Cigar::Match(len) => {
                        let len = *len as usize;
                        for _ in 0..len {
                            let &b = xm_iter.next().unwrap();
                            if pos >= start_1 {
                                break;
                            }
                            if b == b'.' {
                                pos += 1;
                                continue;
                            }

                            let (m, u, c) = match b {
                                b'X' | b'x' => chg.entry(pos).or_insert((0, 0, 0)),
                                b'H' | b'h' => chh.entry(pos).or_insert((0, 0, 0)),
                                b'Z' | b'z' => cpg.entry(pos).or_insert((0, 0, 0)),
                                b'U' | b'u' => unknown.entry(pos).or_insert((0, 0, 0)),
                                _ => return Err(format!("invalid XM key {}", b as char)),
                            };

                            if matches!(b, b'X' | b'H' | b'Z' | b'U') {
                                *m += 1;
                            } else {
                                *u += 1;
                            }
                            *c += 1;
                            pos += 1;
                        }
                    }
                    Cigar::Ins(len) | Cigar::SoftClip(len) => {
                        for _ in 0..*len as usize {
                            xm_iter.next();
                        }
                    }
                    Cigar::Del(len) => {
                        pos += *len as u64;
                    }
                    _ => return Err(format!("invalid CIGAR operation {}", c.char())),
                }
            }
        }

        Ok(())
    }

    pub fn cpg(&self) -> &BTreeMap<String, BTreeMap<u64, (u32, u32, u32)>> {
        &self.cpg
    }

    pub fn chg(&self) -> &BTreeMap<String, BTreeMap<u64, (u32, u32, u32)>> {
        &self.chg
    }

    pub fn chh(&self) -> &BTreeMap<String, BTreeMap<u64, (u32, u32, u32)>> {
        &self.chh
    }

    pub fn chrs(&self) -> &[String] {
        &self.chrs
    }
}

fn end_cigar(start: u64, cigar: &[Cigar]) -> u64 {
    let mut end = start - 1;
    for c in cigar {
        match c {
            Cigar::Match(len) | Cigar::Del(len) | Cigar::RefSkip(len) => end += *len as u64,
            Cigar::Ins(_) | Cigar::SoftClip(_) => {}
            _ => {
                panic!("invalid CIGAR operation: {}", c)
            }
        }
    }
    end
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn process_paired() {
        let chrs = vec![String::from("15")];
        let mut citosyne_genome = CytosineGenome::new(chrs);

        let pos_1 = 59205706;
        let chr_1 = String::from("15");
        let xm_1 = "HhH..HhH..x..h..x..h.h.....xh..x...xh..hh....x..x....Z......x.......x.......h...................x.hh.......hhhh.h.....h.h...........................";
        let cigar_1 = vec![Cigar::Match(7), Cigar::Del(3), Cigar::Match(141)]; // 7M3D141M

        let pos_2 = 59205706;
        let chr_2 = String::from("15");
        let xm_2 = "HhH..HhH..x..h..x..h.h.....xh..x...xh..hh....x..x....Z......x.......x.......h...................x.hh.......hhhh.h.....h.h...........................";
        let cigar_2 = vec![Cigar::Match(7), Cigar::Del(3), Cigar::Match(141)]; // 7M3D141M

        let reverse_2 = false;

        citosyne_genome
            .process_paired(
                pos_1,
                &chr_1,
                xm_1.as_bytes(),
                &cigar_1,
                pos_2,
                &chr_2,
                xm_2.as_bytes(),
                &cigar_2,
                reverse_2,
            )
            .unwrap();

        let pos = 59205762 as u64;
        let (m, u, cov) = citosyne_genome.cpg().get("15").unwrap().get(&pos).unwrap();
        assert_eq!((*m, *u, *cov), (1, 0, 1));
    }

    #[test]
    fn process_paired_2() {
        let chrs = vec![String::from("10")];
        let mut citosyne_genome = CytosineGenome::new(chrs);

        let start_1 = 20706917;
        let chr_1 = String::from("10");
        let xm_1 = "x..x....h...h....h..................h...h.h.hh.......z.....x....z.hh....x.h...................................x....hhh...h.........x..x...h.....Zx..Z".as_bytes();
        let cigar_1 = vec![Cigar::Match(149)]; // 149M

        let start_2 = 20706906;
        let chr_2 = String::from("10");
        let xm_2 = "H......h...x..x....h...h....h..................h...h.h.hh.......z.....x....z.hh....x.h...................................x....hhh...h.........x..x...h".as_bytes();
        let cigar_2 = vec![Cigar::Match(150)]; // 150M

        let reverse_2 = false;

        citosyne_genome
            .process_paired(
                start_1, &chr_1, xm_1, &cigar_1, start_2, &chr_2, xm_2, &cigar_2, reverse_2,
            )
            .unwrap();

        let cpg = citosyne_genome.cpg().get("10").unwrap();

        let pos = 20706970 as u64;
        let (m, u, cov) = cpg.get(&pos).unwrap();
        assert_eq!((*m, *u, *cov), (0, 1, 1));

        let pos = 20706981 as u64;
        let (m, u, cov) = cpg.get(&pos).unwrap();
        assert_eq!((*m, *u, *cov), (0, 1, 1));

        let pos = 20707061 as u64;
        let (m, u, cov) = cpg.get(&pos).unwrap();
        assert_eq!((*m, *u, *cov), (1, 0, 1));

        let pos = 20707065 as u64;
        let (m, u, cov) = cpg.get(&pos).unwrap();
        assert_eq!((*m, *u, *cov), (1, 0, 1));
    }

    #[test]
    fn process_paired_overlap() {
        let chrs = vec![String::from("10")];
        let mut citosyne_genome = CytosineGenome::new(chrs);

        let start_1 = 42092485;
        let chr_1 = String::from("10");
        let xm_1 = ".z...hh........hh.......h...Zx....h...hh....Z.....Z...h....h...h........Z...hh....Z...hh..............hh...hh....h.......z....h....Zx..hh..........Z...".as_bytes();
        let cigar_1 = vec![Cigar::Match(151)]; // 151M

        let start_2 = 42092585;
        let chr_2 = String::from("10");
        let xm_2 = "...uU...HH...hh.......Z....h....Zx..hh..........Z...hh........hh..............Zx........Zx.......Z...hh....Z..h.h.......Z...hh....Z...hh...........h..".as_bytes();
        let cigar_2 = vec![Cigar::Match(2), Cigar::Ins(1), Cigar::Match(147)]; // 2M1I147M

        let reverse_2 = false;

        citosyne_genome
            .process_paired(
                start_1, &chr_1, xm_1, &cigar_1, start_2, &chr_2, xm_2, &cigar_2, reverse_2,
            )
            .unwrap();

        let cpg = citosyne_genome.cpg().get("10").unwrap();

        assert_eq!(cpg.len(), 9);

        let pos = 42092486 as u64;
        let (m, u, cov) = cpg.get(&pos).unwrap();
        assert_eq!((*m, *u, *cov), (0, 1, 1));
        let pos = 42092513 as u64;
        let (m, u, cov) = cpg.get(&pos).unwrap();
        assert_eq!((*m, *u, *cov), (1, 0, 1));
        let pos = 42092529 as u64;
        let (m, u, cov) = cpg.get(&pos).unwrap();
        assert_eq!((*m, *u, *cov), (1, 0, 1));
        let pos = 42092535 as u64;
        let (m, u, cov) = cpg.get(&pos).unwrap();
        assert_eq!((*m, *u, *cov), (1, 0, 1));
        let pos = 42092557 as u64;
        let (m, u, cov) = cpg.get(&pos).unwrap();
        assert_eq!((*m, *u, *cov), (1, 0, 1));
        let pos = 42092567 as u64;
        let (m, u, cov) = cpg.get(&pos).unwrap();
        assert_eq!((*m, *u, *cov), (1, 0, 1));
        let pos = 42092606 as u64;
        let (m, u, cov) = cpg.get(&pos).unwrap();
        assert_eq!((*m, *u, *cov), (0, 1, 1));
        let pos = 42092616 as u64;
        let (m, u, cov) = cpg.get(&pos).unwrap();
        assert_eq!((*m, *u, *cov), (1, 0, 1));
        let pos = 42092632 as u64;
        let (m, u, cov) = cpg.get(&pos).unwrap();
        assert_eq!((*m, *u, *cov), (1, 0, 1));
    }

    #[test]
    fn process_paired_reverse() {
        let chrs = vec![String::from("10")];
        let mut citosyne_genome = CytosineGenome::new(chrs);

        let start_1 = 100110304;
        let chr_1 = String::from("10");
        let xm_1 = "....hh..........h.h.......................h..h...x......x.....h.h.Z.hx......hhx..h...........hh....hx.......h.hx......x.......Z....hx..hx...hh.........".as_bytes();
        let cigar_1 = vec![Cigar::Match(151)]; // 151M

        let start_2 = 100110405;
        let chr_2 = String::from("10");
        let xm_2 = ".......h.hx......x.......Z....hx..hx...hh............hhhh..h.h..h........h.........x....Z.......x...h.hx......hhx..h..h.x......x.....x.........H...HX.".as_bytes();
        let cigar_2 = vec![Cigar::Match(146), Cigar::Del(2), Cigar::Match(4)]; // 146M2D4M

        let reverse_2 = true;

        citosyne_genome
            .process_paired(
                start_1, &chr_1, xm_1, &cigar_1, start_2, &chr_2, xm_2, &cigar_2, reverse_2,
            )
            .unwrap();

        let cpg = citosyne_genome.cpg().get("10").unwrap();

        assert_eq!(cpg.len(), 3);

        let pos = 100110370 as u64;
        let (m, u, cov) = cpg.get(&pos).unwrap();
        assert_eq!((*m, *u, *cov), (1, 0, 1));
        let pos = 100110430 as u64;
        let (m, u, cov) = cpg.get(&pos).unwrap();
        assert_eq!((*m, *u, *cov), (1, 0, 1));
        let pos = 100110493 as u64;
        let (m, u, cov) = cpg.get(&pos).unwrap();
        assert_eq!((*m, *u, *cov), (1, 0, 1));
    }
}

pub fn write_cytosine_report_cpg(
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

pub fn write_cytosine_report_chg(
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
                "{}\t{}\t+\t{}\t{}\tCHG\t{}",
                id,
                pos,
                m,
                u,
                std::str::from_utf8(triplet).unwrap()
            )?;

            let (m, u, _) = chg.get(&(pos + 2)).unwrap_or(&(0, 0, 0));

            writeln!(
                writer,
                "{}\t{}\t-\t{}\t{}\tCHG\t{}",
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

pub fn write_cytosine_report_chh(
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

pub fn write_bed_graph(
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

pub fn write_bismark_cov(
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
