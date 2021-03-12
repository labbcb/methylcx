use std::{collections::BTreeMap, usize};

use rust_htslib::bam::record::Cigar;
struct Context {
    max_len: usize,
    methylated: Vec<u32>,
    unmethylated: Vec<u32>,
    total_methylated: u32,
    total_unmethylated: u32,
}

impl Context {
    fn new(capacity: usize) -> Self {
        Context {
            max_len: 0,
            methylated: Vec::with_capacity(capacity),
            unmethylated: Vec::with_capacity(capacity),
            total_methylated: 0,
            total_unmethylated: 0,
        }
    }

    fn resize(&mut self, new_len: usize) {
        if new_len > self.max_len {
            self.methylated.resize(new_len, 0);
            self.unmethylated.resize(new_len, 0);
            self.max_len = new_len;
        }
    }

    fn add_methylated(&mut self, idx: usize, count: u32) {
        self.methylated[idx] += count;
        self.total_methylated += count;
    }

    fn add_unmethylated(&mut self, idx: usize, count: u32) {
        self.unmethylated[idx] += count;
        self.total_unmethylated += count;
    }

    fn methylated(&self) -> &[u32] {
        &self.methylated
    }

    fn unmethylated(&self) -> &[u32] {
        &self.unmethylated
    }

    fn total_methylated(&self) -> u32 {
        self.total_methylated
    }

    fn total_unmethylated(&self) -> u32 {
        self.total_unmethylated
    }

    fn total(&self) -> u32 {
        self.total_methylated() + self.total_unmethylated()
    }
}

// . for bases not involving cytosines
// X for methylated C in CHG context (was protected)
// x for not methylated C in CHG context (was converted)
// H for methylated C in CHH context (was protected)
// h for not methylated C in CHH context (was converted)
// Z for methylated C in CpG context (was protected)
// z for not methylated C in CpG context (was converted)
// U for methylated C in Unknown context (was protected)
// u for not methylated C in Unknown context (was converted)

pub struct CytosineRead {
    max_len: usize,
    cpg: Context,
    chg: Context,
    chh: Context,
}

impl CytosineRead {
    pub fn new(capacity: usize) -> Self {
        CytosineRead {
            max_len: 0,
            cpg: Context::new(capacity),
            chg: Context::new(capacity),
            chh: Context::new(capacity),
        }
    }

    pub fn process(&mut self, xm: &[u8], reverse: bool, cigar: &[Cigar]) -> Result<(), String> {
        self.resize(xm.len());

        let mut cigar = cigar.to_owned();
        let mut xm = xm.to_owned();
        if reverse {
            xm.reverse();
            cigar.reverse();
        }

        let mut idx = 0;
        let mut xm_iter = xm.iter();
        for c in cigar.into_iter() {
            match c {
                Cigar::Match(len) => {
                    let len = len as usize;
                    for _ in 0..len {
                        let &b = xm_iter.next().unwrap();

                        match b {
                            b'X' => self.chg.add_methylated(idx, 1),
                            b'x' => self.chg.add_unmethylated(idx, 1),
                            b'H' => self.chh.add_methylated(idx, 1),
                            b'h' => self.chh.add_unmethylated(idx, 1),
                            b'Z' => self.cpg.add_methylated(idx, 1),
                            b'z' => self.cpg.add_unmethylated(idx, 1),
                            _ => {}
                        }

                        idx += 1;
                    }
                }
                Cigar::SoftClip(len) => {
                    let len = len as usize;
                    for _ in 0..len {
                        xm_iter.next();
                    }
                }
                Cigar::Ins(len) => {
                    let len = len as usize;
                    for _ in 0..len {
                        xm_iter.next();
                    }
                    idx += len;
                }
                Cigar::Del(_) => {}
                _ => return Err(format!("invalid CIGAR operation {}", c.char())),
            }
        }
        Ok(())
    }

    pub fn resize(&mut self, new_len: usize) {
        if new_len > self.max_len {
            self.cpg.resize(new_len);
            self.chg.resize(new_len);
            self.chh.resize(new_len);
            self.max_len = new_len;
        }
    }

    pub fn cpg_m(&self) -> &[u32] {
        &self.cpg.methylated()
    }

    pub fn cpg_u(&self) -> &[u32] {
        &self.cpg.unmethylated()
    }

    pub fn chg_m(&self) -> &[u32] {
        &self.chg.methylated()
    }

    pub fn chg_u(&self) -> &[u32] {
        &self.chg.unmethylated()
    }

    pub fn chh_m(&self) -> &[u32] {
        &self.chh.methylated()
    }

    pub fn chh_u(&self) -> &[u32] {
        &self.chh.unmethylated()
    }

    pub fn total_cpg_methylated(&self) -> u32 {
        self.cpg.total_methylated()
    }

    pub fn total_cpg_unmethylated(&self) -> u32 {
        self.cpg.total_unmethylated()
    }

    pub fn total_cpg(&self) -> u32 {
        self.cpg.total()
    }

    pub fn total_chg_methylated(&self) -> u32 {
        self.chg.total_methylated()
    }

    pub fn total_chg_unmethylated(&self) -> u32 {
        self.chg.total_unmethylated()
    }

    pub fn total_chg(&self) -> u32 {
        self.chg.total()
    }

    pub fn total_chh_methylated(&self) -> u32 {
        self.chh.total_methylated()
    }

    pub fn total_chh_unmethylated(&self) -> u32 {
        self.chh.total_unmethylated()
    }

    pub fn total_chh(&self) -> u32 {
        self.chh.total()
    }

    pub fn total_c(&self) -> u32 {
        self.total_cpg() + self.total_chg() + self.total_chh()
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

pub struct ClipperSingleConfig {
    pub five_prime_clip: u32,
    pub three_prime_clip: u32,
    pub five_soft_clip: u32,
    pub three_soft_clip: u32,
}

pub struct ClipperPairedConfig {
    pub five_prime_clip_1: u32,
    pub three_prime_clip_1: u32,
    pub five_soft_clip_1: u32,
    pub three_soft_clip_1: u32,
    pub five_prime_clip_2: u32,
    pub three_prime_clip_2: u32,
    pub five_soft_clip_2: u32,
    pub three_soft_clip_2: u32,
}

pub struct ClipperSingle {
    config: ClipperSingleConfig,
    total_m: u32,
    total_sm: u32,
    total_ms: u32,
    total_sms: u32,
}

// pub struct ClipperPaired {
//     config: ClipperPairedConfig,
//     total_m: u32,
//     total_sm: u32,
//     total_ms: u32,
//     total_sms: u32,
// }

/// Clipper trims first and/or last bases from a sequence and updates CIGAR string and POS value.
/// When five_prime_clip is larger than 0 it trims bases from the begining of sequence.
/// WHen three_prime_clip is larger than 0 it trims bases from the end of sequence.
/// If sequence was aligned in reverse mode then values of five_prime_clip and three_prime_clip are swapped.
/// CIGAR operations are considered while trimming bases.
/// The algorithm walks to CIGAR operations until there are bases to clip.
/// Matches (M) consume bases to clip, update POS and trim sequence.
/// Insertions (I) consume bases to clip and trim sequence. POS is not altered.
/// Deletions (D) update POS. Sequence is not trimmed and bases to clip is not clipped.
impl ClipperSingle {
    pub fn new(config: ClipperSingleConfig) -> ClipperSingle {
        ClipperSingle {
            config,
            total_m: 0,
            total_sm: 0,
            total_ms: 0,
            total_sms: 0,
        }
    }

    pub fn process(
        &mut self,
        cigar: Vec<Cigar>,
        pos: u64,
        xm: Vec<u8>,
        reverse: bool,
    ) -> Option<(Vec<Cigar>, u64, Vec<u8>)> {
        self.count_cigar(&cigar);

        let five_prime_clip = if reverse {
            self.config.three_prime_clip
        } else {
            self.config.five_prime_clip
        };
        let three_prime_clip = if reverse {
            self.config.five_prime_clip
        } else {
            self.config.three_prime_clip
        };

        let mut cigar = cigar;
        let mut pos = pos;
        let mut xm = xm;
        if five_prime_clip > 0 {
            match clip_five_prime(cigar, pos, xm, five_prime_clip) {
                Some((new_cigar, new_pos, new_xm)) => {
                    cigar = new_cigar;
                    pos = new_pos;
                    xm = new_xm;
                }
                None => return None,
            }
        }

        if three_prime_clip > 0 {
            match clip_three_prime(cigar, xm, three_prime_clip) {
                Some((new_cigar, new_xm)) => {
                    cigar = new_cigar;
                    xm = new_xm;
                }
                None => return None,
            }
        }

        let five_soft_clip = if reverse {
            self.config.three_soft_clip
        } else {
            self.config.five_soft_clip
        };
        let three_soft_clip = if reverse {
            self.config.five_soft_clip
        } else {
            self.config.three_soft_clip
        };

        if five_soft_clip > 0 {
            match clip_five_soft(cigar, pos, xm, five_soft_clip) {
                Some((new_cigar, new_pos, new_xm)) => {
                    cigar = new_cigar;
                    pos = new_pos;
                    xm = new_xm;
                }
                None => return None,
            }
        }

        if three_soft_clip > 0 {
            match clip_three_soft(cigar, xm, three_soft_clip) {
                Some((new_cigar, new_xm)) => {
                    cigar = new_cigar;
                    xm = new_xm;
                }
                None => return None,
            }
        }

        Some((cigar, pos, xm))
    }

    pub fn total_m(&self) -> u32 {
        self.total_m
    }

    pub fn total_sm(&self) -> u32 {
        self.total_sm
    }

    pub fn total_ms(&self) -> u32 {
        self.total_ms
    }

    pub fn total_sms(&self) -> u32 {
        self.total_sms
    }

    pub fn total(&self) -> u32 {
        self.total_m + self.total_sm + self.total_ms + self.total_sms
    }

    fn count_cigar(&mut self, cigar: &[Cigar]) {
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

        if leading_softclips == 0 && trailing_softclips == 0 {
            self.total_m += 1;
        } else if leading_softclips > 0 && trailing_softclips == 0 {
            self.total_sm += 1;
        } else if leading_softclips == 0 && trailing_softclips > 0 {
            self.total_ms += 1;
        } else {
            self.total_sms += 1;
        }
    }
}

fn clip_five_prime(
    cigar: Vec<Cigar>,
    pos: u64,
    xm: Vec<u8>,
    clip: u32,
) -> Option<(Vec<Cigar>, u64, Vec<u8>)> {
    if let Some(Cigar::Match(_)) = cigar.first() {
        let mut new_cigar: Vec<Cigar> = Vec::new();
        let mut pos = pos;
        let mut start = 0;
        let mut clip = clip;

        let mut cigar_iter = cigar.into_iter();
        while let Some(c) = cigar_iter.next() {
            if clip == 0 {
                new_cigar.push(c);
                continue;
            }
            match c {
                Cigar::Match(len) => {
                    if len > clip {
                        new_cigar.push(Cigar::Match(len - clip));
                        start += clip;
                        pos += clip as u64;
                        clip = 0;
                    } else {
                        start += len;
                        pos += len as u64;
                        clip -= len;
                    }
                }
                Cigar::Ins(len) => {
                    if len > clip {
                        new_cigar.push(Cigar::Ins(len - clip));
                        start += clip;
                        clip = 0;
                    } else {
                        start += len;
                        clip -= len;
                    }
                }
                Cigar::Del(len) => {
                    if len > clip {
                        new_cigar.push(Cigar::Del(len - clip));
                        pos += clip as u64;
                    } else {
                        pos += len as u64;
                    }
                }
                _ => {}
            }
        }

        let start = start as usize;
        if start < xm.len() {
            Some((new_cigar, pos, xm[start..].to_vec()))
        } else {
            None
        }
    } else {
        Some((cigar, pos, xm))
    }
}

fn clip_three_prime(cigar: Vec<Cigar>, xm: Vec<u8>, clip: u32) -> Option<(Vec<Cigar>, Vec<u8>)> {
    if let Some(Cigar::Match(_)) = cigar.last() {
        let mut new_cigar: Vec<Cigar> = Vec::new();
        let mut end = xm.len() as i32;
        let mut clip = clip;

        let mut cigar_iter = cigar.into_iter().rev();
        while let Some(c) = cigar_iter.next() {
            if clip == 0 {
                new_cigar.push(c);
                continue;
            }
            match c {
                Cigar::Match(len) => {
                    if len > clip {
                        new_cigar.push(Cigar::Match(len - clip));
                        end -= clip as i32;
                        clip = 0;
                    } else {
                        end -= len as i32;
                        clip -= len;
                    }
                }
                Cigar::Ins(len) => {
                    if len > clip {
                        new_cigar.push(Cigar::Ins(len - clip));
                        end -= clip as i32;
                        clip = 0;
                    } else {
                        end -= len as i32;
                        clip -= len;
                    }
                }
                Cigar::Del(len) => {
                    if len > clip {
                        new_cigar.push(Cigar::Del(len - clip));
                    }
                }
                _ => {}
            }
        }

        if end > 0 {
            new_cigar.reverse();
            let end = end as usize;
            Some((new_cigar, xm[..end].to_vec()))
        } else {
            None
        }
    } else {
        Some((cigar, xm))
    }
}

fn clip_five_soft(
    cigar: Vec<Cigar>,
    pos: u64,
    xm: Vec<u8>,
    clip: u32,
) -> Option<(Vec<Cigar>, u64, Vec<u8>)> {
    if let Some(Cigar::SoftClip(_)) = cigar.first() {
        let mut new_cigar: Vec<Cigar> = Vec::new();
        let mut pos = pos;
        let mut start = 0;
        let mut clip = clip;

        let mut cigar_iter = cigar.into_iter();
        while let Some(c) = cigar_iter.next() {
            if clip == 0 {
                new_cigar.push(c);
                continue;
            }
            match c {
                Cigar::Match(len) => {
                    if len > clip {
                        new_cigar.push(Cigar::Match(len - clip));
                        start += clip;
                        pos += clip as u64;
                        clip = 0;
                    } else {
                        start += len;
                        pos += len as u64;
                        clip -= len;
                    }
                }
                Cigar::Ins(len) => {
                    if len > clip {
                        new_cigar.push(Cigar::Ins(len - clip));
                        start += clip;
                        clip = 0;
                    } else {
                        start += len;
                        clip -= len;
                    }
                }
                Cigar::Del(len) => {
                    if len > clip {
                        new_cigar.push(Cigar::Del(len - clip));
                        pos += clip as u64;
                    } else {
                        pos += len as u64;
                    }
                }
                _ => {}
            }
        }

        let start = start as usize;
        if start < xm.len() {
            Some((new_cigar, pos, xm[start..].to_vec()))
        } else {
            None
        }
    } else {
        Some((cigar, pos, xm))
    }
}

fn clip_three_soft(cigar: Vec<Cigar>, xm: Vec<u8>, clip: u32) -> Option<(Vec<Cigar>, Vec<u8>)> {
    if let Some(Cigar::SoftClip(_)) = cigar.last() {
        let mut new_cigar: Vec<Cigar> = Vec::new();
        let mut end = xm.len() as i32;
        let mut clip = clip;

        let mut cigar_iter = cigar.into_iter().rev();
        while let Some(c) = cigar_iter.next() {
            if clip == 0 {
                new_cigar.push(c);
                continue;
            }
            match c {
                Cigar::Match(len) => {
                    if len > clip {
                        new_cigar.push(Cigar::Match(len - clip));
                        end -= clip as i32;
                        clip = 0;
                    } else {
                        end -= len as i32;
                        clip -= len;
                    }
                }
                Cigar::Ins(len) => {
                    if len > clip {
                        new_cigar.push(Cigar::Ins(len - clip));
                        end -= clip as i32;
                        clip = 0;
                    } else {
                        end -= len as i32;
                        clip -= len;
                    }
                }
                Cigar::Del(len) => {
                    if len > clip {
                        new_cigar.push(Cigar::Del(len - clip));
                    }
                }
                _ => {}
            }
        }

        if end > 0 {
            new_cigar.reverse();
            let end = end as usize;
            Some((new_cigar, xm[..end].to_vec()))
        } else {
            None
        }
    } else {
        Some((cigar, xm))
    }
}

// impl ClipperPaired {
//     pub fn new(config: ClipperPairedConfig) -> ClipperPaired {
//         ClipperPaired {
//             config,
//             total_m: 0,
//             total_sm: 0,
//             total_ms: 0,
//             total_sms: 0,
//         }
//     }

//     pub fn process(
//         &mut self,
//         cigar_1: Vec<Cigar>,
//         pos_1: u64,
//         xm_1: Vec<u8>,
//         reverse_1: bool,
//         cigar_2: Vec<Cigar>,
//         pos_2: u64,
//         xm_2: Vec<u8>,
//         reverse_2: bool,
//     ) -> Option<((Vec<Cigar>, u64, Vec<u8>), (Vec<Cigar>, u64, Vec<u8>))> {
//         Some(((cigar_1, pos_1, xm_1), (cigar_2, pos_2, xm_2)))
//     }

//     pub fn total_m(&self) -> u32 {
//         self.total_m
//     }

//     pub fn total_sm(&self) -> u32 {
//         self.total_sm
//     }

//     pub fn total_ms(&self) -> u32 {
//         self.total_ms
//     }

//     pub fn total_sms(&self) -> u32 {
//         self.total_sms
//     }

//     pub fn total(&self) -> u32 {
//         self.total_m + self.total_sm + self.total_ms + self.total_sms
//     }
// }

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
