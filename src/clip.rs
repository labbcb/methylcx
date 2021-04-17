use rust_htslib::bam::record::Cigar;

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
