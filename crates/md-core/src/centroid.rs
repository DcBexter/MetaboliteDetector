//! Centroided and nominal-mass spectra, like the old `gcms::SpectrumFunctions::getCentroidSpectrum<int,float>`
//! and `Spectrum<int,float>::getBinnedSpectrum` that the NetCDF export uses (algorithms.md §4.2).

use crate::Chromatogram;
use crate::settings::GcmsSettings;

/// A centroided spectrum: one entry per detected peak.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct Centroids {
    pub mz: Vec<f32>,
    /// Peak area above its baseline, truncated to an integer like the old `int` intensity.
    pub intensity: Vec<i32>,
    /// Number of points at or above `SF_FWHM` times the peak maximum (the old "FWHM" is a point count).
    pub fwhm: Vec<f64>,
}

/// Thresholds of `getCentroidSpectrum`: slope that starts a peak, slope that ends it, and the
/// allowed baseline deviation in standard deviations; `fwhm` (`SF_FWHM`) is the share of the peak
/// maximum a point must reach to count in the centroid mass.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CentroidSettings {
    pub threshold_begin: f64,
    pub threshold_end: f64,
    pub max_base_dev: f64,
    pub fwhm: f64,
}

impl From<&GcmsSettings> for CentroidSettings {
    fn from(s: &GcmsSettings) -> Self {
        Self {
            threshold_begin: s.sf_centroid_threshold_begin,
            threshold_end: s.sf_centroid_threshold_end,
            max_base_dev: s.sf_centroid_max_base_dev,
            fwhm: s.sf_fwhm,
        }
    }
}

impl Default for CentroidSettings {
    fn default() -> Self {
        Self {
            threshold_begin: 10.0,
            threshold_end: -5.0,
            max_base_dev: 30.0,
            fwhm: 0.5,
        }
    }
}

/// Segment length of `analyzeBaseline` for spectra.
const SEGMENT: usize = 25;

/// Old `SavitzkyGolayFilterFirstDerivative<int>`: slope per point, truncated to an integer;
/// the two points at each end are 0.
fn slope(v: &[i32]) -> Vec<i32> {
    let mut d = vec![0; v.len()];
    for (k, w) in v.windows(5).enumerate() {
        let [a, b, c, e, f] = [w[0], w[1], w[2], w[3], w[4]].map(f64::from);
        let _ = c;
        d[k + 2] = ((-2.0 * a - b + e + 2.0 * f) / 10.0) as i32;
    }
    d
}

struct Baseline {
    mean: f64,
    sd: f64,
}

/// Old `SpectrumFunctions::analyzeBaseline<int>`: pooled mean and standard deviation of the
/// non-overlapping segments of 25 points that cross their own mean at least 12 times.
/// A spectrum shorter than a segment makes the old code fall back to segments of one point, whose
/// variance is 0/0: the deviation is NaN, so no peak end is ever accepted. With no accepted
/// segment at all, mean and deviation are 0.
fn analyze_baseline(v: &[i32]) -> Baseline {
    if v.len() < SEGMENT {
        return Baseline {
            mean: 0.0,
            sd: f64::NAN,
        };
    }
    let (mut count, mut mean, mut var) = (0usize, 0.0f64, 0.0f64);
    let mut i = 0;
    while i + SEGMENT < v.len() {
        let seg: Vec<f64> = v[i..i + SEGMENT].iter().map(|x| f64::from(*x)).collect();
        i += SEGMENT;
        let m = seg.iter().sum::<f64>() / SEGMENT as f64;
        let crossings = seg
            .windows(2)
            .filter(|w| (m <= w[0] && w[1] < m) || (w[0] <= m && m < w[1]))
            .count();
        if crossings < SEGMENT / 2 {
            continue;
        }
        let seg_var = seg.iter().map(|x| (x - m).powi(2)).sum::<f64>() / (SEGMENT - 1) as f64;
        if count == 0 {
            (mean, var) = (m, seg_var);
        } else {
            // Pooled mean and sample variance of all accepted points (Chan et al.).
            let (n1, n2) = (count as f64, SEGMENT as f64);
            let pooled = (n1 * mean + n2 * m) / (n1 + n2);
            var = ((n1 - 1.0) * var
                + (n2 - 1.0) * seg_var
                + n1 * (mean - pooled).powi(2)
                + n2 * (m - pooled).powi(2))
                / (n1 + n2 - 1.0);
            mean = pooled;
        }
        count += SEGMENT;
    }
    Baseline {
        mean,
        sd: var.sqrt(),
    }
}

/// Old `SpectrumFunctions::getCentroidSpectrum<int,float>` with the peak search `detectPeaks<int>`:
/// a peak starts where the slope exceeds `threshold_begin`, ends where it stays above
/// `threshold_end` twice, and is widened by its own width on both sides before its baseline,
/// area and centre of mass are taken. Spectra of fewer than 5 points have no peaks.
pub fn centroid(mz: &[f32], intensity: &[i32], s: &CentroidSettings) -> Centroids {
    let n = intensity.len().min(mz.len());
    let (v, mz) = (&intensity[..n], &mz[..n]);
    let mut out = Centroids::default();
    if n < 5 {
        return out;
    }
    let d = slope(v);
    let base = analyze_baseline(v);
    let mut i = 2;
    while i < n - 2 {
        if f64::from(d[i]) > s.threshold_begin {
            let (start, mut top) = (i, i);
            i += 1;
            while i < n - 2 {
                if d[i] >= 0 {
                    top = i;
                    i += 1;
                } else if d[i + 1] < 0 {
                    i += 1;
                    break;
                } else {
                    i += 1;
                }
            }
            let mut end = 0;
            // The search for the end moves the outer position too: the next peak starts after it.
            while i + 1 < n {
                let rising = |k: usize| f64::from(d[k]) > s.threshold_end;
                if rising(i)
                    && rising(i + 1)
                    && (base.sd == 0.0
                        || base.sd * s.max_base_dev > (f64::from(v[i - 1]) - base.mean).abs())
                {
                    end = i - 1;
                    break;
                }
                i += 1;
            }
            if end != 0 && end - start > 2 && top != start && top != end {
                peak(mz, v, start, end, top, s.fwhm, &mut out);
            }
        }
        i += 1;
    }
    out
}

/// `createSpectrumPeak`, `determineIntegral` and `determineCentroid` for the points `first..=last`;
/// `top` is the highest point the slope search found. The old area is the plain sum of the points:
/// the peak baseline is still unset when it is taken.
fn peak(
    mz: &[f32],
    v: &[i32],
    first: usize,
    last: usize,
    top: usize,
    fwhm: f64,
    out: &mut Centroids,
) {
    let w = &v[first..=last];
    let max = (first..=last).fold(top, |m, k| if v[k] > v[m] { k } else { m });
    let least = fwhm * f64::from(v[max]);
    let (mut weighted, mut total, mut count) = (0.0f64, 0.0f64, 0.0);
    for (k, y) in w.iter().enumerate() {
        if f64::from(*y) >= least {
            weighted += f64::from(mz[first + k] * *y as f32);
            total += f64::from(*y);
            count += 1.0;
        }
    }
    out.mz.push((weighted / total) as f32);
    out.intensity
        .push(w.iter().map(|y| f64::from(*y)).sum::<f64>() as i32);
    out.fwhm.push(count);
}

/// Old `Spectrum<int,float>::getBinnedSpectrum`: the highest intensity within each bin of `width`
/// (a bin spans its centre ± `width / 2`), labelled with the bin centre.
pub fn nominal(mz: &[f32], intensity: &[i32], width: f32) -> (Vec<f32>, Vec<i32>) {
    let n = intensity.len().min(mz.len());
    let (mut bins, mut values) = (Vec::new(), Vec::new());
    if n == 0 {
        return (bins, values);
    }
    let half = f64::from(width) / 2.0;
    let centre = |m: f32| (((f64::from(m) + half) / f64::from(width)) as i32) as f32 * width;
    let (mut bin, mut best) = (centre(mz[0]), intensity[0]);
    for i in 0..n {
        if f64::from(mz[i]) >= f64::from(bin) + half {
            bins.push(bin);
            values.push(best);
            best = intensity[i];
            bin = centre(mz[i]);
        } else if intensity[i] > best {
            best = intensity[i];
        }
        if i == n - 1 {
            bins.push(bin);
            values.push(best);
        }
    }
    (bins, values)
}

/// What the NetCDF export writes per scan (the old "export type").
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum ExportKind {
    /// Every point as read.
    #[default]
    Profile,
    /// One point per detected peak.
    Centroid,
    /// The highest intensity within each nominal mass (bins of 1).
    Nominal,
}

impl ExportKind {
    pub const ALL: [ExportKind; 3] = [Self::Profile, Self::Centroid, Self::Nominal];

    pub fn label(self) -> &'static str {
        match self {
            Self::Profile => "Profile",
            Self::Centroid => "Centroid",
            Self::Nominal => "Nominal mass",
        }
    }

    /// Suffix the old export added to the file name.
    pub fn suffix(self) -> &'static str {
        match self {
            Self::Profile => "_profile",
            Self::Centroid => "_centroid",
            Self::Nominal => "_nominal",
        }
    }
}

/// The run as the export writes it: each spectrum centroided or binned. Intensities are cut to
/// integers first, like the old `int` spectra. A scan without peaks stays, with no points.
pub fn reduce(c: &Chromatogram, kind: ExportKind, s: &CentroidSettings) -> Chromatogram {
    if kind == ExportKind::Profile {
        return c.clone();
    }
    let mut out = Chromatogram {
        scan_times_s: c.scan_times_s.clone(),
        offsets: vec![0],
        mz: Vec::new(),
        intensity: Vec::new(),
        exact_masses: c.exact_masses && kind == ExportKind::Centroid,
    };
    for scan in 0..c.scan_count() {
        let (mz, y) = c.spectrum(scan);
        let y: Vec<i32> = y.iter().map(|v| *v as i32).collect();
        match kind {
            ExportKind::Centroid => {
                let k = centroid(mz, &y, s);
                out.mz.extend(k.mz);
                out.intensity.extend(k.intensity.iter().map(|v| *v as f32));
            }
            _ => {
                let (m, v) = nominal(mz, &y, 1.0);
                out.mz.extend(m);
                out.intensity.extend(v.iter().map(|v| *v as f32));
            }
        }
        out.offsets.push(out.mz.len());
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn nominal_takes_the_highest_value_per_bin() {
        let mz = [72.9, 73.2, 73.4, 74.6, 75.0];
        let v = [5, 9, 7, 3, 4];
        let (bins, values) = nominal(&mz, &v, 1.0);
        assert_eq!(bins, [73.0, 75.0]);
        assert_eq!(values, [9, 4]);
        assert!(nominal(&[], &[], 1.0).0.is_empty());
    }

    #[test]
    fn a_gaussian_gives_one_centroid_at_its_centre() {
        let mz: Vec<f32> = (0..60).map(|i| 100.0 + 0.1 * i as f32).collect();
        let v: Vec<i32> = (0..60)
            .map(|i| (1000.0 * (-((i as f64 - 30.0) / 4.0).powi(2) / 2.0).exp()) as i32)
            .collect();
        let c = centroid(&mz, &v, &CentroidSettings::default());
        assert_eq!(c.mz.len(), 1);
        assert!((c.mz[0] - 103.0).abs() < 0.01, "{:?}", c.mz);
        assert!(c.intensity[0] > 5000);
    }

    #[test]
    fn reducing_a_run_keeps_the_scans() {
        let c = Chromatogram {
            scan_times_s: vec![1.0, 2.0],
            offsets: vec![0, 3, 5],
            mz: vec![72.9, 73.2, 74.6, 73.0, 73.4],
            intensity: vec![5.0, 9.0, 3.0, 4.0, 6.0],
            exact_masses: false,
        };
        let n = reduce(&c, ExportKind::Nominal, &CentroidSettings::default());
        assert_eq!(n.offsets, [0, 2, 3]);
        assert_eq!(n.intensity, [9.0, 3.0, 6.0]);
        let k = reduce(&c, ExportKind::Centroid, &CentroidSettings::default());
        assert_eq!((k.scan_count(), k.offsets.clone()), (2, vec![0, 0, 0]));
        assert_eq!(reduce(&c, ExportKind::Profile, &Default::default()), c);
    }

    #[test]
    fn short_spectra_have_no_peaks() {
        let c = centroid(&[1.0, 2.0, 3.0], &[1, 5, 1], &CentroidSettings::default());
        assert!(c.mz.is_empty());
    }
}
