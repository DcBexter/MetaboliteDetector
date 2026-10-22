//! The MID report of a finished scan, without egui: the runs (columns), the substances with their
//! MIDs, the per-run statistics and the text exports. The MID window moves runs and substances and
//! shows the result; a move changes the order only, the values stay as determined.

use crate::batch::pick;
use crate::library::LibraryCompound;

use super::{Method, Mids, group_stat};

/// A move in the result view: run columns (old "Change column order…") or substances with their
/// diagrams (old "Change order of diagrams…"), from one position to another.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Move {
    Run(usize, usize),
    Compound(usize, usize),
}

/// The result of a MID scan: what the tables, plots and exports show.
pub struct MidReport {
    pub method: Method,
    /// The substances with fragments the MIDs refer to.
    pub lib: Vec<LibraryCompound>,
    /// One per fit: the labelled runs' names, or the group.
    pub columns: Vec<String>,
    /// Full path of each per-run column, for its hover; empty for the group methods.
    pub paths: Vec<String>,
    /// Position of each column when determined, before the user moved runs.
    pub runs: Vec<usize>,
    pub mids: Mids,
    /// Per-run methods: mean and standard error of the labelled runs per substance, fragment and
    /// isotopomer (old group `GL`); empty for the group methods.
    pub stats: Vec<Vec<Vec<(f64, f64)>>>,
    /// The old program's text exports: the MIDs, and for the per-run methods the statistics.
    pub text: String,
    pub stats_text: Option<String>,
}

/// The positions `0..n` with `a` and `b` swapped: the order a run move leaves.
fn swapped(n: usize, a: usize, b: usize) -> Vec<usize> {
    let mut order: Vec<usize> = (0..n).collect();
    order.swap(a, b);
    order
}

impl MidReport {
    /// Applies a move to the table, the diagrams and the exports. The values stay as determined:
    /// the statistics keep summing the runs in their first order, so not even the last bit moves.
    pub fn apply(&mut self, m: Move) {
        match m {
            Move::Run(a, b) => {
                let order = swapped(self.columns.len(), a, b);
                self.columns = pick(&self.columns, &order);
                self.runs = pick(&self.runs, &swapped(self.runs.len(), a, b));
                if !self.paths.is_empty() {
                    self.paths = pick(&self.paths, &swapped(self.paths.len(), a, b));
                }
                for f in self
                    .mids
                    .compounds
                    .iter_mut()
                    .flat_map(|c| &mut c.fragments)
                {
                    f.fits = pick(&f.fits, &swapped(f.fits.len(), a, b));
                }
            }
            Move::Compound(a, b) => {
                self.mids.compounds.swap(a, b);
                if !self.stats.is_empty() {
                    self.stats.swap(a, b);
                }
                if self.stats_text.is_some() {
                    let mut first = self.mids.clone();
                    for f in first.compounds.iter_mut().flat_map(|c| &mut c.fragments) {
                        let mut fits = vec![None; f.fits.len()];
                        for (j, fit) in std::mem::take(&mut f.fits).into_iter().enumerate() {
                            fits[self.runs[j]] = fit;
                        }
                        f.fits = fits;
                    }
                    self.stats_text = Some(first.stats_text(&self.lib));
                }
            }
        }
        self.text = self.mids_text();
    }

    /// Statistics and text exports of the MIDs as determined.
    pub fn refresh(&mut self) {
        let per_run = self.method.per_run();
        self.runs = (0..self.columns.len()).collect();
        self.stats = (self.mids.compounds.iter())
            .filter(|_| per_run)
            .map(|c| {
                (c.fragments.iter())
                    .map(|f| (0..f.isotopomers).map(|k| group_stat(f, k)).collect())
                    .collect()
            })
            .collect();
        self.text = self.mids_text();
        self.stats_text = per_run.then(|| self.mids.stats_text(&self.lib));
    }

    pub fn mids_text(&self) -> String {
        match self.method.per_run() {
            true => (self.mids).run_text(&self.lib, self.method, &self.columns),
            false => self.mids.group_text(&self.lib),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::mid::{CompoundMid, FragmentMid, MidFit};

    /// Two substances over the runs `columns`, with the shares `0.1 * (c * 3 + j + 1)` per run.
    fn report(method: Method, columns: &[&str]) -> MidReport {
        let fit = |x: f64| {
            Some(MidFit {
                mids: vec![x, 1.0 - x],
                r2: Some(x),
                cis: vec![Some(x / 10.0), None],
            })
        };
        let compound = |c: usize| CompoundMid {
            compound: c,
            fragments: vec![FragmentMid {
                lo: 73,
                hi: 74,
                isotopomers: 2,
                fits: (0..columns.len())
                    .map(|j| fit(0.1 * (c * 3 + j + 1) as f64))
                    .collect(),
            }],
        };
        let mut r = MidReport {
            method,
            lib: ["Alanine", "Glycine"].map(LibraryCompound::new).to_vec(),
            columns: columns.iter().map(|c| c.to_string()).collect(),
            paths: Vec::new(),
            runs: Vec::new(),
            mids: Mids {
                compounds: vec![compound(0), compound(1)],
                errors: Vec::new(),
            },
            stats: Vec::new(),
            text: String::new(),
            stats_text: None,
        };
        r.refresh();
        r
    }

    /// The statistics as bit patterns, so a change in the last bit shows.
    fn bits(s: &[Vec<Vec<(f64, f64)>>]) -> Vec<Vec<Vec<(u64, u64)>>> {
        (s.iter())
            .map(|c| {
                (c.iter())
                    .map(|f| f.iter().map(|(m, e)| (m.to_bits(), e.to_bits())).collect())
                    .collect()
            })
            .collect()
    }

    #[test]
    fn exports_follow_the_moved_runs_and_substances() {
        let mut d = report(Method::Formula, &["a", "b"]);
        let before = (d.text.clone(), d.stats_text.clone());
        d.apply(Move::Run(1, 0));
        assert_eq!(d.columns, ["b", "a"]);
        let head = "\tb\ta\t\tCI (b)\tCI (a)\t\tR2 (b)\tR2 (a)\n";
        assert!(d.text.starts_with(head), "{}", d.text);
        // Run b of Alanine is 0.2, run a 0.1: the values move with their column.
        assert!(
            d.text
                .contains("73 (M0)\t0.2\t0.1\t\t0.02\t0.01\t\t0.2\t0.1\n")
        );
        d.apply(Move::Compound(0, 1));
        assert!(d.text.starts_with("\tb\ta"));
        assert!(d.text.find("Glycine").unwrap() < d.text.find("Alanine").unwrap());
        let stats = d.stats_text.as_deref().unwrap();
        assert!(stats.find("Glycine").unwrap() < stats.find("Alanine").unwrap());
        // The table's mean of Glycine M0 (0.4 and 0.5) moved up with it.
        assert!((d.stats[0][0][0].0 - 0.45).abs() < 1e-12);
        d.apply(Move::Compound(1, 0));
        d.apply(Move::Run(0, 1));
        assert_eq!((d.text, d.stats_text), before);
    }

    /// A run move leaves the statistics and their text bit for bit. A substance move after it gives
    /// what the same substance move gives without the run move, so the text and the table agree.
    #[test]
    fn moving_runs_keeps_the_statistics_bit_for_bit() {
        // Alanine M0 is 0.1, 0.2, 0.3: the mean is 0.20000000000000004, summed backwards 0.19999999999999998.
        let mut d = report(Method::Formula, &["a", "b", "c"]);
        let (stats, stats_text) = (d.stats.clone(), d.stats_text.clone());
        d.apply(Move::Run(2, 0));
        assert_eq!(d.columns, ["c", "b", "a"]);
        assert_eq!((bits(&d.stats), &d.stats_text), (bits(&stats), &stats_text));
        d.apply(Move::Compound(0, 1));
        let mut plain = report(Method::Formula, &["a", "b", "c"]);
        plain.apply(Move::Compound(0, 1));
        assert_eq!(bits(&d.stats), bits(&plain.stats));
        assert_eq!(d.stats_text, plain.stats_text);
        d.apply(Move::Compound(1, 0));
        d.apply(Move::Run(2, 0));
        assert_eq!(d.columns, ["a", "b", "c"]);
        assert_eq!((bits(&d.stats), d.stats_text), (bits(&stats), stats_text));
    }
}
