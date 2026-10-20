//! Noise, baseline and peak detection per ion trace: the old `DefaultPeakDetector<int,float>`
//! (docs/reverse-engineering/algorithms.md §2–3). Deviations are listed in docs/architecture.md.

use rayon::prelude::*;

use crate::Chromatogram;

/// The old `AN_*` settings; defaults are the old static values (file-formats.md §5.1).
/// `AN_MAX_BASE_DEV` is not read by the old peak detector, so it has no field here.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PeakSettings {
    /// `AN_PEAK_THRESHOLD_BEGIN`: smoothed slope (intensity per scan) that starts a peak.
    pub threshold_begin: f32,
    /// `AN_PEAK_THRESHOLD_END`: the falling flank ends once two slopes in a row exceed this.
    pub threshold_end: f32,
    /// `AN_MIN_PEAK_HEIGHT`: minimum height above the baseline in noise units (S/N).
    pub min_peak_height: f32,
    /// `AN_NUMBER_BINS`: sub-scan steps of the precise apex.
    pub number_bins: u32,
}

impl Default for PeakSettings {
    fn default() -> Self {
        Self {
            threshold_begin: 10.0,
            threshold_end: -5.0,
            min_peak_height: 20.0,
            number_bins: 10,
        }
    }
}

/// A peak in one ion trace. Scan numbers are indices into the run; the intensities are
/// `IonTraces::trace(row)[first_scan..=last_scan]` (for `detect_peaks_matched`: the filtered trace).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Peak {
    /// Row in `IonTraces` (`IonTraces::masses[row]` is the nominal mass).
    pub row: u32,
    pub first_scan: u32,
    /// First scan with the highest intensity.
    pub apex_scan: u32,
    pub last_scan: u32,
    /// Apex of a parabola through `apex_scan ± 1`, in steps of `1 / number_bins` scan.
    /// `None` where the old code refused it (shift of a whole scan or more, apex at the peak edge).
    pub precise_apex: Option<f32>,
    /// Height above the robust baseline divided by the noise units `sqrt(apex) * noise_factor`.
    pub signal_to_noise: f32,
    /// Sum of slopes against the peak shape in percent of those along it (0 = clean peak).
    pub quality: f32,
    pub sharpness: f32,
}

/// Ion-major copy of a run, built once: one intensity trace per nominal mass, so peak detection
/// reads contiguous memory instead of searching every spectrum per ion.
#[derive(Debug, Clone, PartialEq)]
pub struct IonTraces {
    /// Nominal mass (m/z rounded) of each row, ascending.
    pub masses: Vec<u32>,
    scans: usize,
    values: Vec<f32>,
}

// shortcut: m/z is rounded to nominal mass (unit-resolution data, like the example runs);
// high-resolution runs need exact-mass traces plus the old `reorderPeaks` merge (algorithms.md §4.2).
const MASS_LIMIT: usize = 1 << 16;

/// Nominal mass `k` for m/z in `[k - 0.5, k + 0.5)`, the same intervals `IonTraces::new` scans.
fn nominal(mz: f32) -> Option<usize> {
    let m = (mz + 0.5).floor();
    (m >= 0.0 && m < MASS_LIMIT as f32).then_some(m as usize)
}

/// Upper bound for the ion-major index (masses × scans × 4 bytes). Real runs need well under 1 GiB;
/// a corrupt file with huge scan or mass counts must fail instead of exhausting memory.
pub const MAX_TRACE_BYTES: usize = 8 << 30;

#[derive(Debug, Clone, Copy, thiserror::Error, PartialEq, Eq)]
#[error("ion index of {masses} masses x {scans} scans needs more than {} GiB", MAX_TRACE_BYTES >> 30)]
pub struct TooLarge {
    pub masses: usize,
    pub scans: usize,
}

impl IonTraces {
    pub fn new(c: &Chromatogram) -> Result<Self, TooLarge> {
        let scans = c.scan_count();
        let present =
            c.mz.par_chunks(1 << 20)
                .map(|chunk| {
                    let mut p = vec![false; MASS_LIMIT];
                    chunk
                        .iter()
                        .filter_map(|m| nominal(*m))
                        .for_each(|m| p[m] = true);
                    p
                })
                .reduce(
                    || vec![false; MASS_LIMIT],
                    |mut a, b| {
                        a.iter_mut().zip(b).for_each(|(x, y)| *x |= y);
                        a
                    },
                );
        let masses: Vec<u32> = (0..MASS_LIMIT as u32)
            .filter(|m| present[*m as usize])
            .collect();
        let mut row_of = vec![0; MASS_LIMIT];
        masses
            .iter()
            .enumerate()
            .for_each(|(r, m)| row_of[*m as usize] = r);
        let too_large = TooLarge {
            masses: masses.len(),
            scans,
        };
        let cells = masses.len().checked_mul(scans).ok_or(too_large)?;
        if cells > MAX_TRACE_BYTES / size_of::<f32>() {
            return Err(too_large);
        }
        let mut values = vec![0.0f32; cells];
        // Each task owns a block of rows and walks the matching m/z range of every scan once.
        const BLOCK: usize = 8;
        if scans > 0 {
            values
                .par_chunks_mut(scans * BLOCK)
                .zip(masses.par_chunks(BLOCK))
                .for_each(|(rows, block)| {
                    let first_row = row_of[block[0] as usize];
                    let (lo, hi) = (block[0] as f32 - 0.5, block[block.len() - 1] as f32 + 0.5);
                    for s in 0..scans {
                        let (m, i) = c.spectrum(s);
                        let a = m.partition_point(|x| *x < lo);
                        for (mz, int) in m[a..].iter().zip(&i[a..]).take_while(|(x, _)| **x < hi) {
                            let row = row_of[nominal(*mz).unwrap()] - first_row;
                            // Usually one point per nominal mass and scan, so f32 sums stay exact.
                            rows[row * scans + s] += *int;
                        }
                    }
                });
        }
        Ok(Self {
            masses,
            scans,
            values,
        })
    }

    pub fn len(&self) -> usize {
        self.masses.len()
    }

    pub fn is_empty(&self) -> bool {
        self.masses.is_empty()
    }

    pub fn scan_count(&self) -> usize {
        self.scans
    }

    /// Intensity of one nominal mass in every scan.
    pub fn trace(&self, row: usize) -> &[f32] {
        &self.values[row * self.scans..(row + 1) * self.scans]
    }

    /// Peaks of every trace, detected in parallel; sorted by row, then scan.
    pub fn detect_peaks(&self, settings: &PeakSettings) -> Vec<Peak> {
        (0..self.len())
            .into_par_iter()
            .map_init(
                || (Vec::new(), Vec::new()),
                |(smooth, slope), row| {
                    let mut peaks = Vec::new();
                    detect_peaks(self.trace(row), settings, smooth, slope, &mut peaks);
                    peaks.iter_mut().for_each(|p| p.row = row as u32);
                    peaks
                },
            )
            .flat_map_iter(|p| p)
            .collect()
    }
}

/// Savitzky-Golay, window 5, quadratic (old `SavitzkyGolayFilter<int>` and `…FirstDerivative`):
/// smoothed values truncated to integers like the old int filter, and the slope per scan.
/// The two values at each end copy their inner neighbour (smoothing) or are 0 (slope).
fn savitzky_golay(v: &[f32], smooth: &mut Vec<f32>, slope: &mut Vec<f32>) {
    let n = v.len();
    smooth.clear();
    smooth.resize(n, 0.0);
    slope.clear();
    slope.resize(n, 0.0);
    for (k, w) in v.windows(5).enumerate() {
        let [a, b, c, d, e] = [w[0], w[1], w[2], w[3], w[4]].map(f64::from);
        smooth[k + 2] = ((-3.0 * (a + e) + 12.0 * (b + d) + 17.0 * c) / 35.0).trunc() as f32;
        slope[k + 2] = ((-2.0 * a - b + d + 2.0 * e) / 10.0) as f32;
    }
    (smooth[0], smooth[1]) = (smooth[2], smooth[2]);
    (smooth[n - 2], smooth[n - 1]) = (smooth[n - 3], smooth[n - 3]);
}

struct Baseline {
    mean: f64,
    sd: f64,
}

/// Old `analyzeBaseline` (§2.1): pooled mean and standard deviation of the non-overlapping
/// segments of `len` scans whose slopes all stay inside the peak thresholds. The old function
/// also returns a noise factor and MAD, which the detector never reads, so they are not computed.
fn analyze_baseline(smooth: &[f32], slope: &[f32], len: usize, s: &PeakSettings) -> Baseline {
    let (mut count, mut mean, mut m2) = (0.0f64, 0.0f64, 0.0f64);
    let mut start = 0;
    while start + len < smooth.len() {
        let r = start..start + len;
        start += len;
        if !slope[r.clone()]
            .iter()
            .all(|d| *d > s.threshold_end && *d < s.threshold_begin)
        {
            continue;
        }
        let seg = &smooth[r];
        let m = seg.iter().map(|v| *v as f64).sum::<f64>() / len as f64;
        let seg_m2: f64 = seg.iter().map(|v| (*v as f64 - m).powi(2)).sum();
        // Chan's pooled mean and variance of all accepted points.
        let (n, delta) = (len as f64, m - mean);
        let total = count + n;
        mean += delta * n / total;
        m2 += seg_m2 + delta * delta * count * n / total;
        count = total;
    }
    Baseline {
        mean,
        sd: if count > 1.0 {
            (m2 / (count - 1.0)).sqrt()
        } else {
            0.0
        },
    }
}

/// Old `DefaultPeakDetector::detectPeaks` (§3.1) on one trace; appends to `out`.
pub fn detect_peaks(
    v: &[f32],
    s: &PeakSettings,
    smooth: &mut Vec<f32>,
    slope: &mut Vec<f32>,
    out: &mut Vec<Peak>,
) {
    let n = v.len();
    if n < 5 {
        return;
    }
    savitzky_golay(v, smooth, slope);
    let (sm, d) = (&smooth[..], &slope[..]);
    let base = analyze_baseline(sm, d, 15, s);

    let mut tailing = false;
    let mut i = 2;
    while i < n - 2 {
        if let Some((beg, en, mx)) = candidate(&mut i, d, sm, &base, s, tailing) {
            match check_peak(v, beg, en, mx, s) {
                Some(p) => {
                    let last = v[p.last_scan as usize] as f64 - base.mean;
                    tailing = last > 0.2 * (v[p.apex_scan as usize] as f64 - base.mean);
                    out.push(p);
                    // The old loop steps back one scan here, so a tailing peak starts where the last one ended + 1.
                    i -= 1;
                }
                None => {
                    tailing = false;
                    i = beg;
                }
            }
        } else {
            tailing = false;
        }
        i += 1;
    }
}

/// Rising and falling flank search from scan `i`; leaves `i` where the old loop left it.
fn candidate(
    i: &mut usize,
    d: &[f32],
    sm: &[f32],
    base: &Baseline,
    s: &PeakSettings,
    tailing: bool,
) -> Option<(usize, usize, usize)> {
    let n = d.len();
    if !(d[*i] > s.threshold_begin || tailing) {
        return None;
    }
    let beg = *i;
    let mut mx = Some(beg);
    let mut flat = 0;
    *i += 1;
    while *i + 1 < n {
        if d[*i] >= 0.0 {
            if d[*i] == 0.0 {
                flat += 1;
                if flat > 2 && (sm[*i] as f64 - base.mean).abs() <= base.sd {
                    mx = None;
                    break;
                }
            }
            mx = Some(*i);
        } else if d[*i + 1] < 0.0 {
            break;
        }
        *i += 1;
    }
    let mx = mx?;
    while *i + 1 < n {
        if d[*i] > s.threshold_end && d[*i + 1] > s.threshold_end {
            let en = *i - 1;
            return (en - beg > 3 && mx != beg && mx != en).then_some((beg, en, mx));
        }
        *i += 1;
    }
    None
}

/// Noise factor of the old detector. The old `GCMSAnalyzer` sets it from the run's global noise factor,
/// which is 1.0 on the example runs (algorithms.md §2.3, docs/architecture.md).
pub const NOISE_FACTOR: f64 = 1.0;

/// Index of the highest value, starting from `start`: a later index wins only with a strictly
/// higher value (strict `>`, like the old code).
fn argmax(w: &[f32], start: usize) -> usize {
    (0..w.len()).fold(start, |m, k| if w[k] > w[m] { k } else { m })
}

/// Old `createPeak` (§3.2, border extension off as in the old default) and `calculateBaseline`
/// (§2.2): fit a robust baseline and keep the peak if it stands `min_peak_height` noise units above it.
fn check_peak(v: &[f32], first: usize, last: usize, mx: usize, s: &PeakSettings) -> Option<Peak> {
    let w = &v[first..=last];
    // createPeak starts the apex search at the last rising scan of the flank search.
    let apex = argmax(w, mx - first);
    let max = w[apex] as f64;
    let at = |k: usize| w[k] as f64;

    // Line through the lowest points left of and from the apex on (both start at the apex, strict
    // `<`), then pivot it around the lower one (right one on a tie) until no point lies more than
    // 0.1 below it. Same arithmetic as the old code: residual ties decide which points the fit uses.
    let lowest = |r: std::ops::Range<usize>| r.fold(apex, |m, k| if w[k] < w[m] { k } else { m });
    let (l, r) = (lowest(0..apex), lowest(apex..w.len()));
    let pivot = if w[l] < w[r] { l } else { r };
    let mut slope = if l == r {
        0.0
    } else {
        (at(r) - at(l)) / (r - l) as f64
    };
    let mut icpt = at(l) - l as f64 * slope;
    for _ in 0..=w.len() {
        let mut moved = false;
        for k in 0..w.len() {
            if (k as f64 * slope + icpt) - at(k) > 0.1 {
                slope = if k < pivot {
                    (at(pivot) - at(k)) / (pivot - k) as f64
                } else {
                    (at(k) - at(pivot)) / (k - pivot) as f64
                };
                icpt = at(pivot) - pivot as f64 * slope;
                moved = true;
            }
        }
        if !moved {
            break;
        }
    }
    // Least-squares line through the first half (+1) of the distinct residuals, lowest first
    // (old std::map<double, index>: a repeated residual keeps its first index and counts once).
    let mut res: Vec<(f64, usize)> = (0..w.len())
        .map(|k| (at(k) - (k as f64 * slope + icpt), k))
        .collect();
    res.sort_by(|a, b| a.0.total_cmp(&b.0).then(a.1.cmp(&b.1)));
    res.dedup_by(|a, b| a.0 == b.0);
    res.truncate(res.len() / 2 + 1);
    let (icpt, fit_slope) = fit_line(res.iter().map(|(_, k)| (*k as f64, at(*k))));

    let nu = max.sqrt() * NOISE_FACTOR;
    let height = max - (icpt + fit_slope * apex as f64);
    if !(nu > 0.0 && height >= s.min_peak_height as f64 * nu) {
        return None;
    }
    Some(peak(
        w,
        first,
        apex,
        NOISE_FACTOR,
        (height / nu) as f32,
        s.number_bins,
    ))
}

/// Old `Peak` constructor (§3.5): a peak over `w = trace[first..first + w.len()]` with its metrics.
fn peak(w: &[f32], first: usize, apex: usize, nf: f64, signal_to_noise: f32, bins: u32) -> Peak {
    Peak {
        row: 0,
        first_scan: first as u32,
        apex_scan: (first + apex) as u32,
        last_scan: (first + w.len() - 1) as u32,
        precise_apex: precise_apex(w, apex, bins).map(|t| t + first as f32),
        signal_to_noise,
        quality: quality(w, apex),
        sharpness: sharpness(w, apex, (w[apex] as f64).sqrt() * nf),
    }
}

/// Old `tools::Filter::matched2ndDerivativeGaussianFilter` (§3.3): correlation with a negative
/// second derivative of a Gaussian, sigma = `density / 4`, cut at ±`density`, zeros outside the trace.
/// It responds most to Gaussian peaks with sigma ≈ `density / 4` scans.
pub fn matched_filter<T: Copy + Into<f64>>(v: &[T], density: usize) -> Vec<f64> {
    let d = density.max(1) as isize;
    let kernel: Vec<f64> = (-d..=d)
        .map(|x| {
            let u = 4.0 * x as f64 / d as f64;
            1.73 / (d as f64).sqrt() * (1.0 - u * u) * (-0.5 * u * u).exp()
        })
        .collect();
    (0..v.len() as isize)
        .map(|i| {
            (-d..=d)
                .filter_map(|x| {
                    v.get(usize::try_from(i + x).ok()?)
                        .map(|y| (x, (*y).into()))
                })
                .map(|(x, y)| kernel[(x + d) as usize] * y)
                .sum()
        })
        .collect()
}

/// Old `MatchedPeakDetector<int,float>::detectPeaks` (§3.3–3.4): regions where the filtered trace
/// (truncated to integers like the old int filter) rises above `threshold`, more than 2 scans long.
/// As in the old code, the peak metrics use the filtered values with noise factor 0.5, there is
/// no baseline check, and `signal_to_noise` stays 0. The old app takes `density` and `threshold`
/// from the RI calibration dialog; it has no defaults.
pub fn detect_peaks_matched(v: &[f32], density: usize, threshold: i32, bins: u32) -> Vec<Peak> {
    let f: Vec<f32> = matched_filter(v, density)
        .iter()
        .map(|x| x.trunc() as f32)
        .collect();
    let thr = threshold as f32;
    let mut peaks = Vec::new();
    let mut i = 1;
    while i < f.len() {
        if f[i - 1] < thr && thr < f[i] {
            let Some(end) = (i + 1..f.len()).find(|j| f[*j] <= thr) else {
                break;
            };
            if end - i > 2 {
                let w = &f[i..=end];
                peaks.push(peak(w, i, argmax(w, 0), 0.5, 0.0, bins));
            }
            i = end;
        }
        i += 1;
    }
    peaks
}

/// Ordinary least squares `y = icpt + slope * x`; a single point gives a flat line.
fn fit_line(points: impl Iterator<Item = (f64, f64)> + Clone) -> (f64, f64) {
    let (n, sx, sy) = points.clone().fold((0.0, 0.0, 0.0), |(n, sx, sy), (x, y)| {
        (n + 1.0, sx + x, sy + y)
    });
    let (mx, my) = (sx / n, sy / n);
    let (sxx, sxy) = points.fold((0.0, 0.0), |(a, b), (x, y)| {
        (a + (x - mx) * (x - mx), b + (x - mx) * (y - my))
    });
    let slope = if sxx > 0.0 { sxy / sxx } else { 0.0 };
    (my - slope * mx, slope)
}

/// §3.6: vertex of the parabola through apex − 1, apex, apex + 1 of the peak values, rounded to
/// `1 / bins` scan, relative to the peak start; `None` at the peak edge, on a flat top or for a shift of
/// ≥ 1 scan. Exact half-step ties round down: the old SVD fit in absolute bins lands just below the tie in
/// 217 of 267 ties of the example runs (docs/architecture.md, "Golden test").
fn precise_apex(v: &[f32], apex: usize, bins: u32) -> Option<f32> {
    if apex == 0 || apex + 1 >= v.len() {
        return None;
    }
    let (a, b, c) = (v[apex - 1] as f64, v[apex] as f64, v[apex + 1] as f64);
    let shift = 0.5 * (a - c) / (a - 2.0 * b + c);
    let bins = bins as f64;
    let t = ((apex as f64 + shift) * bins - 0.5).ceil() / bins;
    ((t - apex as f64).abs() < 1.0).then_some(t as f32)
}

/// §3.7: slopes against the peak shape in percent of slopes along it.
fn quality(w: &[f32], apex: usize) -> f32 {
    let (mut along, mut against) = (0.0f64, 0.0f64);
    for (k, p) in w.windows(2).enumerate() {
        let step = (p[1] - p[0]) as f64;
        let rising = if k < apex { step } else { -step };
        if rising >= 0.0 {
            along += rising;
        } else {
            against -= rising;
        }
    }
    if along == 0.0 {
        100.0
    } else {
        (100.0 * against / along) as f32
    }
}

/// §3.8: steepest drop from the apex per scan and noise unit, averaged over both sides.
/// Each side starts at -1 like the old code, so a peak without one side scores lower.
fn sharpness(w: &[f32], apex: usize, nu: f64) -> f32 {
    let max = w[apex] as f64;
    let side = |ks: &mut dyn Iterator<Item = usize>| {
        ks.filter(|k| w[*k] >= 0.0)
            .map(|k| (max - w[k] as f64) / (k.abs_diff(apex) as f64 * nu))
            .fold(-1.0, f64::max)
    };
    let s = (side(&mut (0..apex)) + side(&mut (apex + 1..w.len()))) / 2.0;
    if s.is_finite() {
        s.max(0.0) as f32
    } else {
        0.0
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn precise_apex_rounds_ties_down() {
        // Vertex at 1.25 scans = 12.5 bins: the old fit mostly lands just below the tie.
        assert_eq!(precise_apex(&[1.0, 4.0, 3.0], 1, 10), Some(1.2));
        assert_eq!(precise_apex(&[3.0, 4.0, 1.0], 1, 10), Some(0.7));
        assert_eq!(precise_apex(&[5.0, 5.0, 5.0], 1, 10), None);
    }

    /// Baseline of 1000 with deterministic noise of about ±20 (sd ≈ 12).
    fn noisy(n: usize) -> Vec<f32> {
        let mut x = 12345u32;
        (0..n)
            .map(|_| {
                x = x.wrapping_mul(1_103_515_245).wrapping_add(12345);
                1000.0 + ((x >> 16) % 41) as f32 - 20.0
            })
            .collect()
    }

    fn gauss(v: &mut [f32], center: f64, sigma: f64, height: f64) {
        for (k, y) in v.iter_mut().enumerate() {
            *y += (height * (-0.5 * ((k as f64 - center) / sigma).powi(2)).exp()) as f32;
        }
    }

    fn detect(v: &[f32]) -> Vec<Peak> {
        let mut out = Vec::new();
        detect_peaks(
            v,
            &PeakSettings::default(),
            &mut Vec::new(),
            &mut Vec::new(),
            &mut out,
        );
        out
    }

    #[test]
    fn finds_separate_gaussian_peaks() {
        let mut v = noisy(400);
        gauss(&mut v, 100.0, 2.5, 50_000.0);
        gauss(&mut v, 250.3, 3.0, 8_000.0);
        let p = detect(&v);
        assert_eq!(p.len(), 2, "{p:?}");
        assert_eq!((p[0].apex_scan, p[1].apex_scan), (100, 250));
        assert!(p[0].first_scan < 95 && p[0].last_scan > 105, "{:?}", p[0]);
        assert!(p[0].signal_to_noise > p[1].signal_to_noise && p[1].signal_to_noise > 20.0);
        let apex = p[1].precise_apex.unwrap();
        assert!((apex - 250.3).abs() <= 0.2, "{apex}");
        assert!(p[0].quality < 10.0 && p[0].sharpness > 0.0, "{:?}", p[0]);
    }

    #[test]
    fn ignores_noise_short_and_empty_traces() {
        assert!(detect(&[]).is_empty());
        assert!(detect(&[0.0, 9.0, 99.0, 9.0]).is_empty());
        assert!(detect(&noisy(500)).is_empty());
        assert!(detect(&vec![0.0; 100]).is_empty());
        let mut small = noisy(300);
        gauss(&mut small, 150.0, 2.5, 60.0); // far below 20 noise units
        assert!(detect(&small).is_empty());
    }

    #[test]
    fn flat_topped_peak_counts_once() {
        let mut v = noisy(300);
        gauss(&mut v, 150.0, 3.0, 40_000.0);
        v[148..=152].fill(41_000.0); // detector saturation
        let p = detect(&v);
        assert_eq!(p.len(), 1, "{p:?}");
        // Like the old code: the apex search starts at the last scan whose slope is still ≥ 0, the
        // plateau's middle, and earlier equal values do not replace it.
        assert_eq!(p[0].apex_scan, 150);
    }

    #[test]
    fn peaks_cut_by_the_trace_edges_do_not_panic() {
        for center in [0.0, 1.0, 3.0, 196.0, 198.0, 199.0] {
            let mut v = noisy(200);
            gauss(&mut v, center, 3.0, 30_000.0);
            for p in detect(&v) {
                assert!(p.first_scan <= p.apex_scan && p.apex_scan <= p.last_scan);
                assert!((p.last_scan as usize) < v.len());
            }
        }
    }

    #[test]
    fn matched_filter_finds_peaks_of_its_width() {
        let mut impulse = [0.0; 9];
        impulse[4] = 1.0;
        let k = matched_filter(&impulse, 8);
        assert!(
            k[4] > 0.0 && k[0] < 0.0 && (k[0] - k[8]).abs() < 1e-12,
            "{k:?}"
        );
        let mut v = noisy(300);
        gauss(&mut v, 80.0, 2.0, 5_000.0);
        gauss(&mut v, 200.0, 2.0, 5_000.0);
        let p = detect_peaks_matched(&v, 8, 2_000, 10);
        assert_eq!(p.len(), 2, "{p:?}");
        assert_eq!((p[0].apex_scan, p[1].apex_scan), (80, 200));
        assert!(detect_peaks_matched(&[], 8, 10, 10).is_empty());
    }

    #[test]
    fn ion_traces_take_odd_mz_from_broken_files() {
        let c = Chromatogram {
            scan_times_s: vec![0.0],
            offsets: vec![0, 4],
            mz: vec![-0.5, -0.4, 0.5, f32::INFINITY],
            intensity: vec![1.0, 2.0, 4.0, 8.0],
        };
        let t = IonTraces::new(&c).unwrap();
        assert_eq!(t.masses, vec![0, 1]);
        assert_eq!((t.trace(0), t.trace(1)), (&[3.0][..], &[4.0][..]));
    }

    #[test]
    fn ion_traces_group_by_nominal_mass_and_detect_in_parallel() {
        let mut big = noisy(200);
        gauss(&mut big, 80.0, 2.5, 30_000.0);
        let scans = big.len();
        let c = Chromatogram {
            scan_times_s: (0..scans).map(|s| s as f64).collect(),
            offsets: (0..=scans).map(|s| 3 * s).collect(),
            mz: (0..scans).flat_map(|_| [72.9f32, 73.2, 147.0]).collect(),
            intensity: (0..scans).flat_map(|s| [1.0, big[s] - 1.0, 5.0]).collect(),
        };
        let t = IonTraces::new(&c).unwrap();
        assert_eq!(t.masses, vec![73, 147]);
        assert_eq!(t.trace(0), &big[..]);
        assert_eq!(t.trace(1), &vec![5.0; scans][..]);
        let p = t.detect_peaks(&PeakSettings::default());
        assert_eq!(p.len(), 1, "{p:?}");
        assert_eq!((p[0].row, p[0].apex_scan), (0, 80));
        let empty = Chromatogram {
            scan_times_s: vec![],
            offsets: vec![0],
            mz: vec![],
            intensity: vec![],
        };
        assert!(
            IonTraces::new(&empty)
                .unwrap()
                .detect_peaks(&PeakSettings::default())
                .is_empty()
        );
    }

    #[test]
    fn ion_traces_refuse_absurd_sizes() {
        // 65 536 masses in the first scan, 40 000 scans: 10 GB of f32 if it were allocated.
        let scans = 40_000;
        let mut offsets = vec![0, 65_536];
        offsets.resize(scans + 1, 65_536);
        let c = Chromatogram {
            scan_times_s: vec![0.0; scans],
            offsets,
            mz: (0..65_536).map(|m| m as f32).collect(),
            intensity: vec![1.0; 65_536],
        };
        assert_eq!(
            IonTraces::new(&c),
            Err(TooLarge {
                masses: 65_536,
                scans
            })
        );
    }
}
