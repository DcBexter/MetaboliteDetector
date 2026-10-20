//! Detail of one component, the old component panel (`Compound::getAveragePeakQuality`,
//! `getAverageSignalToNoise`, mean `getSharpness`) and the comparison of its spectrum with a library
//! spectrum (normal / split / difference), see docs/reverse-engineering/algorithms.md §7g.

use std::io;

use crate::csv::write_columns;
use crate::deconv::Component;
use crate::library::LibrarySpectrum;
use crate::peaks::Peak;
use crate::spectrum::{Normalization, normalize};

/// One ion peak of a component.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PeakRow {
    pub mass: u32,
    pub quality: f32,
    pub sharpness: f32,
    pub signal_to_noise: f32,
    pub integral: f64,
}

/// Values of the component panel.
#[derive(Debug, Clone, PartialEq)]
pub struct Detail {
    /// One row per ion mass, sharpest first: the old `map<float, PeakEssential*>` keeps the first of
    /// two peaks with the same mass (§4.5).
    pub peaks: Vec<PeakRow>,
    /// Mean S/N of `peaks`; 0 without peaks.
    pub mean_signal_to_noise: f64,
    /// Mean sharpness of `peaks` (the old tooltip says "highest", the code takes the mean); 0 without peaks.
    pub mean_sharpness: f64,
    /// "Peak dis.indx": mean peak quality, 0 is best; -1 without peaks like the old code.
    pub discrepancy: f64,
}

/// `all[p]` is peak `p` of the list the components were deconvoluted from, `masses[p]` its ion mass.
pub fn component_detail(c: &Component, all: &[Peak], masses: &[u32]) -> Detail {
    let mut peaks: Vec<PeakRow> = Vec::new();
    for f in &c.peaks {
        // Components read from an old `.cmp` come without a peak list: their ions have no statistics.
        let (Some(p), Some(&mass)) = (all.get(f.peak), masses.get(f.peak)) else {
            continue;
        };
        if peaks.iter().all(|r| r.mass != mass) {
            peaks.push(PeakRow {
                mass,
                quality: p.quality,
                sharpness: p.sharpness,
                signal_to_noise: p.signal_to_noise,
                integral: f.integral,
            });
        }
    }
    let mean = |f: fn(&PeakRow) -> f32| match peaks.len() {
        0 => None,
        n => Some(peaks.iter().map(|r| f64::from(f(r))).sum::<f64>() / n as f64),
    };
    Detail {
        mean_signal_to_noise: mean(|r| r.signal_to_noise).unwrap_or(0.0),
        mean_sharpness: mean(|r| r.sharpness).unwrap_or(0.0),
        discrepancy: mean(|r| r.quality).unwrap_or(-1.0),
        peaks,
    }
}

/// The found and the library spectrum on a common mass axis, both normalised on their own.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct Comparison {
    /// Nominal masses, ascending; every mass of either spectrum.
    pub mz: Vec<f64>,
    pub found: Vec<f64>,
    pub library: Vec<f64>,
}

/// Aligns `found` (nominal mass, height) with the library spectrum (m/z rounded to a nominal mass, equal
/// masses added up) and normalises each by `norm`. Without normalisation the base peak is used: the found
/// heights are detector counts and the library is scaled to 0..999, so raw values cannot be compared.
pub fn compare(found: &[(u32, u32)], lib: &LibrarySpectrum, norm: Normalization) -> Comparison {
    let norm = match norm {
        Normalization::None => Normalization::BasePeak,
        n => n,
    };
    let mut f: Vec<(i64, f64)> = found.iter().map(|&(m, h)| (m.into(), h.into())).collect();
    let mut l: Vec<(i64, f64)> = (lib.masses.iter().zip(&lib.values))
        .map(|(m, v)| (m.round() as i64, f64::from(*v)))
        .collect();
    let mut mz: Vec<i64> = f.iter().chain(&l).map(|p| p.0).collect();
    mz.sort_unstable();
    mz.dedup();
    let place = |s: &mut Vec<(i64, f64)>| -> Vec<f64> {
        s.sort_by_key(|p| p.0);
        let mut out = vec![0.0; mz.len()];
        for (m, v) in s.iter() {
            out[mz.binary_search(m).unwrap()] += v;
        }
        normalize(&out, norm)
    };
    Comparison {
        found: place(&mut f),
        library: place(&mut l),
        mz: mz.into_iter().map(|m| m as f64).collect(),
    }
}

impl Comparison {
    /// Found minus library, the third view of the old split plot.
    pub fn difference(&self) -> Vec<f64> {
        self.found
            .iter()
            .zip(&self.library)
            .map(|(f, l)| f - l)
            .collect()
    }

    /// Both spectra and their difference as CSV (file and clipboard).
    pub fn write_csv(&self, mut w: impl io::Write) -> io::Result<()> {
        let diff = self.difference();
        let head = ["m/z", "Found", "Library", "Difference"];
        write_columns(
            &mut w,
            &head,
            &[&self.mz, &self.found, &self.library, &diff],
        )
    }

    pub fn csv_text(&self) -> String {
        let mut out = Vec::new();
        self.write_csv(&mut out).expect("writing to memory");
        String::from_utf8(out).expect("CSV is ASCII")
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::deconv::FittedPeak;

    fn peak(sn: f32, quality: f32, sharpness: f32) -> Peak {
        Peak {
            row: 0,
            first_scan: 0,
            apex_scan: 1,
            last_scan: 2,
            precise_apex: None,
            signal_to_noise: sn,
            quality,
            sharpness,
        }
    }

    fn component(fitted: &[(usize, f64)]) -> Component {
        Component {
            first_scan: 0,
            apex_scan: 1,
            last_scan: 2,
            model: Vec::new(),
            spectrum: Vec::new(),
            peaks: fitted
                .iter()
                .map(|&(peak, integral)| FittedPeak {
                    peak,
                    model_share: 1.0,
                    integral,
                })
                .collect(),
            exact_spectrum: Vec::new(),
        }
    }

    #[test]
    fn means_run_over_one_peak_per_mass() {
        let all = [
            peak(10.0, 2.0, 8.0),
            peak(99.0, 99.0, 99.0),
            peak(20.0, 4.0, 4.0),
        ];
        // Peaks 0 and 1 are both ion 73: the first (sharpest) one counts.
        let d = component_detail(
            &component(&[(0, 5.0), (1, 6.0), (2, 7.0)]),
            &all,
            &[73, 73, 116],
        );
        assert_eq!(d.peaks.len(), 2);
        assert_eq!((d.peaks[0].mass, d.peaks[0].integral), (73, 5.0));
        assert_eq!(d.mean_signal_to_noise, 15.0);
        assert_eq!(d.mean_sharpness, 6.0);
        assert_eq!(d.discrepancy, 3.0);
    }

    #[test]
    fn peaks_missing_from_the_list_are_skipped() {
        let d = component_detail(&component(&[(0, 5.0)]), &[], &[73]);
        assert!(d.peaks.is_empty() && d.discrepancy == -1.0);
    }

    #[test]
    fn a_component_without_peaks_keeps_the_old_defaults() {
        let d = component_detail(&component(&[]), &[], &[]);
        assert_eq!(
            (d.mean_signal_to_noise, d.mean_sharpness, d.discrepancy),
            (0.0, 0.0, -1.0)
        );
    }

    fn library(masses: &[f32], values: &[i32]) -> LibrarySpectrum {
        let mut l = crate::library::LibraryCompound::new("x").spectrum;
        (l.masses, l.values) = (masses.to_vec(), values.to_vec());
        l
    }

    #[test]
    fn spectra_share_one_mass_axis_and_differ_by_subtraction() {
        let c = compare(
            &[(73, 50), (147, 100)],
            &library(&[73.2, 116.0], &[10, 20]),
            Normalization::BasePeak,
        );
        assert_eq!(c.mz, [73.0, 116.0, 147.0]);
        assert_eq!(c.found, [0.5, 0.0, 1.0]);
        assert_eq!(c.library, [0.5, 1.0, 0.0]);
        assert_eq!(c.difference(), [0.0, -1.0, 1.0]);
    }

    #[test]
    fn csv_has_both_spectra_and_the_difference() {
        // "None" compares on the base peak: counts and a 0..999 library are not on one scale.
        let lib = library(&[73.0], &[999]);
        let c = compare(&[(73, 4000), (147, 2000)], &lib, Normalization::None);
        assert_eq!(
            c.csv_text().lines().collect::<Vec<_>>(),
            ["m/z,Found,Library,Difference", "73,1,1,0", "147,0.5,0,0.5"]
        );
    }
}
