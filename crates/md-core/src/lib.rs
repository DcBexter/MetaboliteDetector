//! Core of the MetaboliteDetector rewrite: GC/MS data model, file import and analysis.
//! No GUI dependencies, so the same code drives the desktop app and a future batch CLI.

use rayon::prelude::*;

pub mod andi;
pub mod batch;
pub mod cmp;
pub mod csv;
pub mod deconv;
pub mod diskscan;
pub mod isotope;
pub mod library;
pub mod library_add;
pub mod normalize;
pub mod peaks;
pub mod results;
pub mod ri;
pub mod search;
pub mod settings;
pub mod spectrum;
pub mod stats;

/// Reads a run by extension: the old MetaboliteDetector cache (`.bin` + `.idx`), else ANDI-MS NetCDF.
pub fn read_run(path: &std::path::Path) -> Result<Chromatogram, String> {
    match path
        .extension()
        .is_some_and(|e| e.eq_ignore_ascii_case("bin"))
    {
        true => diskscan::read(path).map_err(|e| e.to_string()),
        false => andi::read(path).map_err(|e| e.to_string()),
    }
}

/// A GC/MS run: one mass spectrum per scan, stored as flat arrays (ANDI-MS layout).
///
/// Values are `f32` like the original (`GCMSDiskScan<int, float>`): half the memory of `f64`
/// for runs with tens of millions of points. Sums are accumulated in `f64`.
#[derive(Debug, Clone, PartialEq)]
pub struct Chromatogram {
    /// Retention time of each scan in seconds, ascending.
    pub scan_times_s: Vec<f64>,
    /// Scan `i` occupies `mz[offsets[i]..offsets[i + 1]]`; `offsets.len() == scan_count() + 1`.
    pub offsets: Vec<usize>,
    /// m/z values, ascending within each scan.
    pub mz: Vec<f32>,
    pub intensity: Vec<f32>,
}

impl Chromatogram {
    pub fn scan_count(&self) -> usize {
        self.scan_times_s.len()
    }

    /// Mass spectrum of one scan as (m/z, intensity) slices.
    pub fn spectrum(&self, scan: usize) -> (&[f32], &[f32]) {
        let r = self.offsets[scan]..self.offsets[scan + 1];
        (&self.mz[r.clone()], &self.intensity[r])
    }

    /// Total ion chromatogram, summed from the spectra (ANDI `total_intensity` is optional).
    pub fn tic(&self) -> Vec<f64> {
        (0..self.scan_count())
            .into_par_iter()
            .map(|s| self.spectrum(s).1.iter().map(|i| *i as f64).sum())
            .collect()
    }

    /// Selected ion chromatogram: per scan, summed intensity of ions within `mz ± tolerance`.
    pub fn sic(&self, mz: f32, tolerance: f32) -> Vec<f64> {
        (0..self.scan_count())
            .into_par_iter()
            .map(|s| {
                let (m, i) = self.spectrum(s);
                let lo = m.partition_point(|x| *x < mz - tolerance);
                let hi = m.partition_point(|x| *x <= mz + tolerance);
                i[lo..hi].iter().map(|v| *v as f64).sum()
            })
            .collect()
    }

    /// Index of the scan whose retention time is closest to `time_s`.
    pub fn scan_at(&self, time_s: f64) -> Option<usize> {
        let idx = self.scan_times_s.partition_point(|t| *t < time_s);
        [idx.checked_sub(1), Some(idx)]
            .into_iter()
            .flatten()
            .filter(|i| *i < self.scan_count())
            .min_by(|a, b| {
                (self.scan_times_s[*a] - time_s)
                    .abs()
                    .total_cmp(&(self.scan_times_s[*b] - time_s).abs())
            })
    }
}

/// Old `SavitzkyGolayFilterFirstDerivative<double>` (window 5, quadratic) as plotted by the old
/// `ChromatogramPlot::setPlotDerivative`: slope per scan, first and last two values 0; runs of
/// fewer than 5 scans come back unchanged, as in the old code.
pub fn first_derivative(v: &[f64]) -> Vec<f64> {
    if v.len() < 5 {
        return v.to_vec();
    }
    let mut out = vec![0.0; v.len()];
    for (k, w) in v.windows(5).enumerate() {
        out[k + 2] = (-2.0 * w[0] - w[1] + w[3] + 2.0 * w[4]) / 10.0;
    }
    out
}

/// Time axis unit of the old "Time scale" setting.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum TimeUnit {
    Milliseconds,
    Seconds,
    #[default]
    Minutes,
}

impl TimeUnit {
    pub const ALL: [TimeUnit; 3] = [Self::Milliseconds, Self::Seconds, Self::Minutes];

    pub fn label(self) -> &'static str {
        match self {
            Self::Milliseconds => "ms",
            Self::Seconds => "s",
            Self::Minutes => "min",
        }
    }

    /// Converts seconds (the run's native unit) into this unit.
    pub fn from_seconds(self, s: f64) -> f64 {
        match self {
            Self::Milliseconds => s * 1000.0,
            Self::Seconds => s,
            Self::Minutes => s / 60.0,
        }
    }

    pub fn to_seconds(self, t: f64) -> f64 {
        match self {
            Self::Milliseconds => t / 1000.0,
            Self::Seconds => t,
            Self::Minutes => t * 60.0,
        }
    }
}

/// Min/max decimation of an x-sorted series to at most `2 * buckets` points within `[x_min, x_max]`.
/// Keeps every peak visible while plotting cost stays proportional to screen width, not data size.
pub fn decimate(xs: &[f64], ys: &[f64], x_min: f64, x_max: f64, buckets: usize) -> Vec<[f64; 2]> {
    let lo = xs.partition_point(|x| *x < x_min).saturating_sub(1);
    let hi = (xs.partition_point(|x| *x <= x_max) + 1).min(xs.len());
    let (xs, ys) = (&xs[lo..hi], &ys[lo..hi]);
    if xs.len() <= 2 * buckets.max(1) {
        return xs.iter().zip(ys).map(|(x, y)| [*x, *y]).collect();
    }
    let per = xs.len().div_ceil(buckets);
    let mut out = Vec::with_capacity(2 * buckets);
    for (cx, cy) in xs.chunks(per).zip(ys.chunks(per)) {
        let (imin, imax) = (0..cy.len()).fold((0, 0), |(a, b), i| {
            (
                if cy[i] < cy[a] { i } else { a },
                if cy[i] > cy[b] { i } else { b },
            )
        });
        let (a, b) = (imin.min(imax), imin.max(imax));
        out.push([cx[a], cy[a]]);
        if b != a {
            out.push([cx[b], cy[b]]);
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn derivative_of_a_ramp_is_its_slope_with_zero_ends() {
        let v: Vec<f64> = (0..8).map(|i| 3.0 * i as f64).collect();
        assert_eq!(
            first_derivative(&v),
            [0.0, 0.0, 3.0, 3.0, 3.0, 3.0, 0.0, 0.0]
        );
        assert_eq!(first_derivative(&[1.0, 2.0, 4.0]), [1.0, 2.0, 4.0]);
    }

    #[test]
    fn time_units_round_trip() {
        for u in TimeUnit::ALL {
            assert!((u.to_seconds(u.from_seconds(90.0)) - 90.0).abs() < 1e-9);
        }
        assert_eq!(TimeUnit::Minutes.from_seconds(90.0), 1.5);
    }

    #[test]
    fn tic_sic_and_scan_lookup() {
        let c = Chromatogram {
            scan_times_s: vec![1.0, 2.0, 3.0],
            offsets: vec![0, 2, 3, 5],
            mz: vec![73.0, 147.0, 73.0, 73.1, 200.0],
            intensity: vec![10.0, 5.0, 20.0, 1.0, 4.0],
        };
        assert_eq!(c.tic(), vec![15.0, 20.0, 5.0]);
        assert_eq!(c.sic(73.0, 0.5), vec![10.0, 20.0, 1.0]);
        assert_eq!(c.spectrum(2), (&[73.1, 200.0][..], &[1.0, 4.0][..]));
        assert_eq!(c.scan_at(2.4), Some(1));
        assert_eq!(c.scan_at(99.0), Some(2));
        assert_eq!(c.scan_at(-5.0), Some(0));
    }

    #[test]
    fn decimate_keeps_extremes_and_bounds_size() {
        let xs: Vec<f64> = (0..10_000).map(f64::from).collect();
        let mut ys = vec![1.0; 10_000];
        ys[4321] = 99.0;
        ys[7000] = -5.0;
        let d = decimate(&xs, &ys, 0.0, 10_000.0, 100);
        assert!(d.len() <= 200);
        assert!(d.contains(&[4321.0, 99.0]) && d.contains(&[7000.0, -5.0]));
        assert!(d.windows(2).all(|w| w[0][0] < w[1][0]));
        // Small inputs and zoomed windows pass through unchanged.
        assert_eq!(decimate(&xs, &ys, 10.0, 12.0, 100).len(), 5);
    }
}
