use rust_htslib::bam::record::Cigar;
use std::io::Write;
use std::{fs::File, io};
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

pub fn write_mbias_single(cytosine_read: &CytosineRead, writer: &mut File) -> io::Result<()> {
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

pub fn write_mbias_paired(
    cytosine_read_1: &CytosineRead,
    cytosine_read_2: &CytosineRead,
    writer: &mut File,
) -> io::Result<()> {
    writer.write_all(b"Strand,Context,Cycle,Methylated,Unmethylated,Coverage\n")?;
    for (i, (m, um)) in cytosine_read_1
        .cpg_m()
        .iter()
        .zip(cytosine_read_1.cpg_u())
        .enumerate()
    {
        writeln!(writer, "Forward,CpG,{},{},{},{}", i + 1, m, um, m + um)?;
    }
    for (i, (m, um)) in cytosine_read_1
        .chg_m()
        .iter()
        .zip(cytosine_read_1.chg_u())
        .enumerate()
    {
        writeln!(writer, "Forward,CHG,{},{},{},{}", i + 1, m, um, m + um)?;
    }
    for (i, (m, um)) in cytosine_read_1
        .chh_m()
        .iter()
        .zip(cytosine_read_1.chh_u())
        .enumerate()
    {
        writeln!(writer, "Forward,CHH,{},{},{},{}", i + 1, m, um, m + um)?;
    }

    for (i, (m, um)) in cytosine_read_2
        .cpg_m()
        .iter()
        .zip(cytosine_read_2.cpg_u())
        .enumerate()
    {
        writeln!(writer, "Reverse,CpG,{},{},{},{}", i + 1, m, um, m + um)?;
    }
    for (i, (m, um)) in cytosine_read_2
        .chg_m()
        .iter()
        .zip(cytosine_read_2.chg_u())
        .enumerate()
    {
        writeln!(writer, "Reverse,CHG,{},{},{},{}", i + 1, m, um, m + um)?;
    }
    for (i, (m, um)) in cytosine_read_2
        .chh_m()
        .iter()
        .zip(cytosine_read_2.chh_u())
        .enumerate()
    {
        writeln!(writer, "Reverse,CHH,{},{},{},{}", i + 1, m, um, m + um)?;
    }

    Ok(())
}
