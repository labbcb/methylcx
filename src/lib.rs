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

    pub fn process(&mut self, xm: &[u8], reverse: bool, cigar: &[Cigar]) {
        let mut cigar = cigar.to_owned();
        let mut xm = xm.to_owned();
        if reverse {
            xm.reverse();
            cigar.reverse();
        }

        let filter: Vec<_> = cigar
            .iter()
            .flat_map(|c| {
                let keep = c.char() != 'S' || c.char() != 'I';
                vec![keep; c.len() as usize]
            })
            .collect();

        let xm: Vec<_> = xm
            .iter()
            .zip(filter)
            .filter(|(_, y)| *y)
            .map(|(x, _)| *x)
            .collect();

        self.resize(xm.len());

        for (i, &b) in xm.iter().enumerate() {
            match b {
                b'X' => self.chg.add_methylated(i, 1),
                b'x' => self.chg.add_unmethylated(i, 1),
                b'H' => self.chh.add_methylated(i, 1),
                b'h' => self.chh.add_unmethylated(i, 1),
                b'Z' => self.cpg.add_methylated(i, 1),
                b'z' => self.cpg.add_unmethylated(i, 1),
                _ => {}
            }
        }
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

    pub fn process(
        &mut self,
        pos: u64,
        chr: String,
        xm: &[u8],
        cigar: &[Cigar],
    ) -> Result<(), String> {
        if !self.chrs.contains(&chr) {
            return Err(format!("invalid chromosome {}", chr));
        }

        let cpg = self.cpg.get_mut(&chr).unwrap();
        let chg = self.chg.get_mut(&chr).unwrap();
        let chh = self.chh.get_mut(&chr).unwrap();
        let unknown = self.unknown.get_mut(&chr).unwrap();

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
                    let len = *len as usize;
                    for _ in 0..len {
                        xm_iter.next();
                    }
                }
                Cigar::Del(len) => {
                    let len = *len as u64;
                    pos += len;
                }
                _ => return Err(format!("invalid CIGAR operation {}", c.char())),
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
}
