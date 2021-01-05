use std::collections::BTreeMap;

use rust_htslib::bam::record::Cigar;

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
    cpg_m: Vec<u32>,
    cpg_u: Vec<u32>,
    chg_m: Vec<u32>,
    chg_u: Vec<u32>,
    chh_m: Vec<u32>,
    chh_u: Vec<u32>,
    total_cpg_m: u32,
    total_cpg_u: u32,
    total_chg_m: u32,
    total_chg_u: u32,
    total_chh_m: u32,
    total_chh_u: u32,
}

impl CytosineRead {
    pub fn new(capacity: usize) -> Self {
        CytosineRead {
            max_len: 0,
            cpg_m: Vec::with_capacity(capacity),
            cpg_u: Vec::with_capacity(capacity),
            chg_m: Vec::with_capacity(capacity),
            chg_u: Vec::with_capacity(capacity),
            chh_m: Vec::with_capacity(capacity),
            chh_u: Vec::with_capacity(capacity),
            total_cpg_m: 0,
            total_cpg_u: 0,
            total_chg_m: 0,
            total_chg_u: 0,
            total_chh_m: 0,
            total_chh_u: 0,
        }
    }

    pub fn process(&mut self, xm: &[u8], reverse: bool, cigar: &[Cigar]) {
        let mut cigar = cigar.to_owned();
        let mut xm = xm.to_owned();
        if reverse {
            xm.reverse();
            cigar.reverse();
        }

        let filter: Vec<bool> = cigar
            .iter()
            .flat_map(|c| {
                let keep = c.char() != 'S' || c.char() != 'I';
                vec![keep; c.len() as usize]
            })
            .collect();

        let xm: Vec<u8> = xm
            .iter()
            .zip(filter)
            .filter(|(_, y)| *y)
            .map(|(x, _)| *x)
            .collect();

        self.resize(xm.len());

        for (i, &b) in xm.iter().enumerate() {
            if b == b'.' {
            } else if b == b'X' {
                self.chg_m[i] += 1;
                self.total_chg_m += 1;
            } else if b == b'x' {
                self.chg_u[i] += 1;
                self.total_chg_u += 1;
            } else if b == b'H' {
                self.chh_m[i] += 1;
                self.total_chh_m += 1;
            } else if b == b'h' {
                self.chh_u[i] += 1;
                self.total_chh_u += 1;
            } else if b == b'Z' {
                self.cpg_m[i] += 1;
                self.total_cpg_m += 1;
            } else if b == b'z' {
                self.cpg_u[i] += 1;
                self.total_cpg_u += 1;
            }
        }
    }

    pub fn resize(&mut self, new_len: usize) {
        if new_len > self.max_len {
            self.cpg_m.resize(new_len, 0);
            self.cpg_u.resize(new_len, 0);
            self.chg_m.resize(new_len, 0);
            self.chg_u.resize(new_len, 0);
            self.chh_m.resize(new_len, 0);
            self.chh_u.resize(new_len, 0);
            self.max_len = new_len;
        }
    }

    pub fn cpg_m(&self) -> &[u32] {
        &self.cpg_m
    }

    pub fn cpg_u(&self) -> &[u32] {
        &self.cpg_u
    }

    pub fn chg_m(&self) -> &[u32] {
        &self.chg_m
    }

    pub fn chg_u(&self) -> &[u32] {
        &self.chg_u
    }

    pub fn chh_m(&self) -> &[u32] {
        &self.chh_m
    }

    pub fn chh_u(&self) -> &[u32] {
        &self.chh_u
    }

    pub fn total_cpg_methylated(&self) -> u32 {
        self.total_cpg_m
    }

    pub fn total_cpg_unmethylated(&self) -> u32 {
        self.total_cpg_u
    }

    pub fn total_cpg(&self) -> u32 {
        self.total_cpg_m + self.total_cpg_u
    }

    pub fn total_chg_methylated(&self) -> u32 {
        self.total_chg_m
    }

    pub fn total_chg_unmethylated(&self) -> u32 {
        self.total_chg_u
    }

    pub fn total_chg(&self) -> u32 {
        self.total_chg_m + self.total_cpg_u
    }

    pub fn total_chh_methylated(&self) -> u32 {
        self.total_chh_m
    }

    pub fn total_chh_unmethylated(&self) -> u32 {
        self.total_chh_u
    }

    pub fn total_chh(&self) -> u32 {
        self.total_chh_m + self.total_chh_u
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

        // Substract 1 base from start position.
        let mut pos = pos - 1;

        let mut xm_iter = xm.iter();
        for c in cigar.iter() {
            match c {
                Cigar::Match(len) => {
                    let len = *len as usize;
                    for _ in 0..len {
                        pos += 1;
                        
                        let &b = xm_iter.next().unwrap();
                        if b == b'.' {
                            continue;
                        }

                        let (m, u, c) = if b == b'X' || b == b'x' {
                            chg.entry(pos).or_insert((0, 0, 0))
                        } else if b == b'H' || b == b'h' {
                            chh.entry(pos).or_insert((0, 0, 0))
                        } else if b == b'Z' || b == b'z' {
                            cpg.entry(pos).or_insert((0, 0, 0))
                        } else if b == b'U' || b == b'u' {
                            unknown.entry(pos).or_insert((0, 0, 0))
                        } else {
                            return Err(format!("invalid XM key {}", b as char));
                        };

                        // Do not count methylation at unknown context.
                        if b == b'X' || b == b'H' || b == b'Z' {
                            *m += 1;
                        } else {
                            *u += 1;
                        }
                        *c += 1;
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
