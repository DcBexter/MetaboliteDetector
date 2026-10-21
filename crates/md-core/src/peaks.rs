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
    /// Peak detection method; `Matched` ignores the three thresholds above.
    pub detector: Detector,
    /// `SF_FWHM`: share of a peak's maximum a point must reach to count in its centroid
    /// ([`IonTraces::new`] takes it; kept here so a change re-runs the analysis).
    pub fwhm: f64,
}

/// "Peak detection method" of the old Calibration Setter window (algorithms.md §3.3). The old main
/// window always uses `Default`; the method is not part of the old settings file.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum Detector {
    /// Old `DefaultPeakDetector`: slope thresholds, robust baseline and S/N (§3.1).
    #[default]
    Default,
    /// Old `MatchedPeakDetector(width, threshold)` ("Matched Filtering"): regions where the trace,
    /// filtered with a matched filter of density `4 * width` scans, lies above `threshold`.
    Matched { width: u32, threshold: i32 },
}

impl Detector {
    /// The old window's values after switching to "Matched Filtering": width 10, threshold 3.
    pub const MATCHED: Self = Self::Matched {
        width: 10,
        threshold: 3,
    };
}

impl Default for PeakSettings {
    fn default() -> Self {
        Self {
            threshold_begin: 10.0,
            threshold_end: -5.0,
            min_peak_height: 20.0,
            number_bins: 10,
            detector: Detector::Default,
            fwhm: 0.5,
        }
    }
}

/// A peak in one ion trace. Scan numbers are indices into the run; the intensities are
/// `IonTraces::trace(row)[first_scan..=last_scan]` (filtered values after a matched detection).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Peak {
    /// Row in `IonTraces` (`IonTraces::masses[row]` is the nominal mass, `IonTraces::mz[row]` the exact one).
    pub row: u32,
    pub first_scan: u32,
    /// First scan with the highest intensity.
    pub apex_scan: u32,
    pub last_scan: u32,
    /// Apex of a parabola through `apex_scan ± 1`, in steps of `1 / number_bins` scan.
    /// `None` where the old code refused it (shift of a whole scan or more, apex at the peak edge);
    /// see [`Peak::interpolated_max_scan`] for what the old code then uses.
    pub precise_apex: Option<f32>,
    /// Height above the robust baseline divided by the noise units `sqrt(apex) * noise_factor`.
    pub signal_to_noise: f32,
    /// Sum of slopes against the peak shape in percent of those along it (0 = clean peak).
    pub quality: f32,
    pub sharpness: f32,
}

impl Peak {
    /// The old `interpolatedMaxScan` field: the precise apex, or 0 where the old code refused it (the
    /// field keeps its initial 0, so deconvolution and `reorderPeaks` put such a peak into bin 0).
    pub fn interpolated_max_scan(&self) -> f64 {
        self.precise_apex.map_or(0.0, f64::from)
    }
}

/// Ion-major copy of a run, built once: one intensity trace per mass, so peak detection reads
/// contiguous memory instead of searching every spectrum per ion.
///
/// `new` makes one trace per nominal mass, like the old CDF import, which centroids every scan to
/// nominal masses. `exact` makes one trace per distinct m/z, like the old `GCMSDiskScan` for runs
/// that keep their exact masses (FastFlight); algorithms.md §4.1.
#[derive(Debug, Clone, PartialEq)]
pub struct IonTraces {
    /// Nominal mass (m/z rounded) of each row, ascending. Repeats in a high-resolution run.
    pub masses: Vec<u32>,
    /// Exact mass of each row, ascending: the distinct m/z of a high-resolution run, else the
    /// nominal mass.
    pub mz: Vec<f32>,
    /// One trace per exact m/z (`IonTraces::exact`).
    exact: bool,
    scans: usize,
    values: Vec<f32>,
}

const MASS_LIMIT: usize = 1 << 16;

/// Nominal mass `k` for m/z in `[k - 0.5, k + 0.5)`, the same intervals `IonTraces::new` scans.
fn nominal(mz: f32) -> Option<usize> {
    let m = (mz + 0.5).floor();
    (m >= 0.0 && m < MASS_LIMIT as f32).then_some(m as usize)
}

/// Distinct m/z of the run within the nominal mass range, ascending.
fn exact_masses(c: &Chromatogram) -> Vec<f32> {
    let mut mz: Vec<f32> =
        c.mz.par_chunks(1 << 20)
            .flat_map_iter(|chunk| {
                let mut v: Vec<f32> = chunk
                    .iter()
                    .copied()
                    .filter(|m| nominal(*m).is_some())
                    .collect();
                v.sort_unstable_by(f32::total_cmp);
                v.dedup();
                v
            })
            .collect();
    mz.par_sort_unstable_by(f32::total_cmp);
    mz.dedup();
    mz
}

/// Nominal masses present in the run, ascending.
fn nominal_masses(c: &Chromatogram) -> Vec<u32> {
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
    (0..MASS_LIMIT as u32)
        .filter(|m| present[*m as usize])
        .collect()
}

/// Gap above which the old centroiding starts a new mass peak.
const CENTROID_GAP: f64 = 0.1;

/// Last point of the centroid part starting at `start`: before an m/z gap above 0.1, or the bottom of
/// a valley after two falls in a row (two rises follow). Past the end the last value read repeats.
fn centroid_part_end(m: &[f32], v: &[f32], start: usize) -> usize {
    let n = m.len();
    let at = |i: usize, missing: f32| if i < n { v[i] } else { missing };
    let mut falling = false;
    let mut i = start;
    while i < n {
        let (v0, m1) = (v[i], if i + 1 < n { m[i + 1] } else { m[i] });
        let v1 = at(i + 1, v0);
        let v2 = at(i + 2, v1);
        if (m1 - m[i]) as f64 > CENTROID_GAP {
            break;
        }
        if !falling && v1 < v0 && v2 < v1 {
            falling = true;
            continue;
        }
        if falling && v1 > v0 && v2 > v1 {
            break;
        }
        i += 1;
    }
    i.min(n - 1)
}

/// The old `Spectrum::getCentroidSpectrum` of one scan (algorithms.md §4.1): points split at m/z
/// gaps above 0.1 and at valleys (two falls, then two rises); each part gives the sum of its points
/// of at least `fwhm` (`SF_FWHM`) times its maximum, at their intensity-weighted m/z.
fn centroid_scan(m: &[f32], v: &[f32], fwhm: f64, out: &mut Vec<(f32, f32)>) {
    let mut i = 0;
    while i < m.len() {
        let start = i;
        i = centroid_part_end(m, v, start);
        let max = v[start..=i].iter().fold(f64::MIN, |a, x| a.max(*x as f64));
        let (mz, sum) = if start == i {
            (m[start] as f64, v[start] as f64)
        } else {
            let (mut num, mut sum) = (0.0f64, 0.0f64);
            for (x, y) in m[start..=i].iter().zip(&v[start..=i]) {
                if *y as f64 >= fwhm * max {
                    sum += *y as f64;
                    num += (*y * *x) as f64;
                }
            }
            let mz = if num == 0.0 && sum == 0.0 {
                -1.0
            } else {
                num / sum
            };
            (mz, sum)
        };
        if mz >= 0.0 {
            out.push((mz as f32, (sum + 0.5).trunc() as f32)); // C `(int)` truncates
        }
        i += 1;
    }
}

/// The run as the old CDF import stores it: every scan centroided. `None` when no scan has two
/// points within 0.1 m/z, where centroiding keeps every point of non-negative integer intensity and m/z >= 0.
/// shortcut: that skip keeps fractional intensities and m/z in [-0.5, 0), which centroiding rounds
/// or drops (CDF intensities are integers); round both ways alike if such input matters.
fn centroided(c: &Chromatogram, fwhm: f64) -> Option<Chromatogram> {
    let close = |s: usize| {
        let (m, _) = c.spectrum(s);
        m.windows(2).any(|p| (p[1] - p[0]) as f64 <= CENTROID_GAP)
    };
    if !(0..c.scan_count()).into_par_iter().any(close) {
        return None;
    }
    let scans: Vec<Vec<(f32, f32)>> = (0..c.scan_count())
        .into_par_iter()
        .map(|s| {
            let (m, i) = c.spectrum(s);
            let mut out = Vec::with_capacity(m.len());
            centroid_scan(m, i, fwhm, &mut out);
            out
        })
        .collect();
    let mut offsets = Vec::with_capacity(scans.len() + 1);
    offsets.push(0);
    for s in &scans {
        offsets.push(offsets.last().unwrap() + s.len());
    }
    let (mz, intensity) = scans.into_iter().flatten().unzip();
    Some(Chromatogram {
        scan_times_s: c.scan_times_s.clone(),
        offsets,
        mz,
        intensity,
        exact_masses: false,
        baseline: c.baseline.clone(),
    })
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
    /// One trace per nominal mass of the centroided scans (`centroided`); centroids of one nominal
    /// mass in a scan are summed. A run with exact masses (FastFlight) gets [`Self::exact`].
    /// `fwhm` is `SF_FWHM` ([`PeakSettings::fwhm`]).
    pub fn new(c: &Chromatogram, fwhm: f64) -> Result<Self, TooLarge> {
        if c.exact_masses {
            return Self::exact(c);
        }
        match centroided(c, fwhm) {
            Some(cc) => Self::build(&cc, false),
            None => Self::build(c, false),
        }
    }

    /// One trace per distinct m/z (the old `getSignal` takes the point of exactly that mass; a
    /// missing mass reads [`Chromatogram::missing_value`]). Masses outside the nominal mass range are dropped.
    pub fn exact(c: &Chromatogram) -> Result<Self, TooLarge> {
        Self::build(c, true)
    }

    fn build(c: &Chromatogram, exact: bool) -> Result<Self, TooLarge> {
        let scans = c.scan_count();
        let (masses, mz): (Vec<u32>, Vec<f32>) = if exact {
            let mz = exact_masses(c);
            (mz.iter().map(|m| nominal(*m).unwrap() as u32).collect(), mz)
        } else {
            let masses = nominal_masses(c);
            let mz = masses.iter().map(|m| *m as f32).collect();
            (masses, mz)
        };
        // Nominal mass -> row (unit resolution only).
        let mut row_of = vec![0; MASS_LIMIT];
        if !exact {
            masses
                .iter()
                .enumerate()
                .for_each(|(r, m)| row_of[*m as usize] = r);
        }
        let too_large = TooLarge {
            masses: masses.len(),
            scans,
        };
        let cells = masses.len().checked_mul(scans).ok_or(too_large)?;
        if cells > MAX_TRACE_BYTES / size_of::<f32>() {
            return Err(too_large);
        }
        let mut values = vec![0.0f32; cells];
        // A mass missing from a scan reads the scan's baseline mean unless it is corrected (§4.1).
        if !c.baseline.corrected && c.baseline.means.len() == scans && scans > 0 {
            (values.par_chunks_mut(scans)).for_each(|row| row.copy_from_slice(&c.baseline.means));
        }
        // Each task owns a block of rows and walks the matching m/z range of every scan once.
        const BLOCK: usize = 8;
        if scans > 0 {
            values
                .par_chunks_mut(scans * BLOCK)
                .zip(mz.par_chunks(BLOCK))
                .for_each(|(rows, block)| {
                    let (first, last) = (block[0], block[block.len() - 1]);
                    for s in 0..scans {
                        let (m, i) = c.spectrum(s);
                        // The first point of a cell replaces the missing value, later ones add up.
                        let mut prev = usize::MAX;
                        let mut put = |cell: usize, v: f32| {
                            match cell == prev {
                                true => rows[cell] += v,
                                false => rows[cell] = v,
                            }
                            prev = cell;
                        };
                        if exact {
                            let a = m.partition_point(|x| *x < first);
                            for (mz, int) in
                                m[a..].iter().zip(&i[a..]).take_while(|(x, _)| **x <= last)
                            {
                                put(block.partition_point(|x| x < mz) * scans + s, *int);
                            }
                            continue;
                        }
                        let first_row = row_of[first as usize];
                        let a = m.partition_point(|x| *x < first - 0.5);
                        for (mz, int) in m[a..]
                            .iter()
                            .zip(&i[a..])
                            .take_while(|(x, _)| **x < last + 0.5)
                        {
                            let row = row_of[nominal(*mz).unwrap()] - first_row;
                            // Usually one point per nominal mass and scan, so f32 sums stay exact.
                            put(row * scans + s, *int);
                        }
                    }
                });
        }
        Ok(Self {
            masses,
            mz,
            exact,
            scans,
            values,
        })
    }

    /// Old `GCMSAnalyzer::detectAllPeaks(true, 10, 0.05)` of the main window and the batch
    /// (algorithms.md §4.0): [`Self::detect_peaks`], then on exact traces [`reorder_peaks`].
    pub fn detect_all_peaks(&mut self, settings: &PeakSettings) -> Vec<Peak> {
        let mut peaks = self.detect_peaks(settings);
        if self.exact {
            reorder_peaks(&mut peaks, self, 10.0, 0.05, settings.number_bins);
        }
        peaks
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

    /// Intensity of one mass (row) in every scan.
    pub fn trace(&self, row: usize) -> &[f32] {
        &self.values[row * self.scans..(row + 1) * self.scans]
    }

    /// Peaks of every trace, detected in parallel; sorted by row, then scan.
    /// With `Detector::Matched` every trace is first replaced by its filtered values (truncated to
    /// integers): the old peaks hold the filtered values, and the deconvolution reads them (§3.3).
    /// So a matched detection is one-shot: a second call would filter the filtered traces.
    pub fn detect_peaks(&mut self, settings: &PeakSettings) -> Vec<Peak> {
        let Detector::Matched { width, threshold } = settings.detector else {
            return self.detect_default(settings);
        };
        if self.scans == 0 {
            return Vec::new();
        }
        let kernel = matched_kernel(4 * width as usize);
        self.values
            .par_chunks_mut(self.scans)
            .enumerate()
            .map_init(Vec::new, |raw, (row, trace)| {
                raw.clear();
                raw.extend_from_slice(trace);
                let mut peaks = Vec::new();
                matched_peaks(
                    raw,
                    &kernel,
                    threshold,
                    settings.number_bins,
                    trace,
                    &mut peaks,
                );
                peaks.iter_mut().for_each(|p| p.row = row as u32);
                peaks
            })
            .flat_map_iter(|p| p)
            .collect()
    }

    fn detect_default(&self, settings: &PeakSettings) -> Vec<Peak> {
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
        precise_apex: precise_apex(w, first, apex, bins),
        signal_to_noise,
        quality: quality(w, apex),
        sharpness: sharpness(w, apex, (w[apex] as f64).sqrt() * nf),
    }
}

/// Old `GCMSAnalyzer::reorderPeaks` (algorithms.md §4.2): peaks on neighbouring exact masses at most
/// `ion_dist` apart whose apexes lie within `scan_dist` scans form one profile, and every member
/// covering the apex of the summed profile takes that apex. Sorts `peaks` by row, then first scan.
/// The old code runs it only for runs not imported from CDF (`DataType != 2`), so
/// [`IonTraces::detect_all_peaks`] calls it only on exact traces.
pub fn reorder_peaks(
    peaks: &mut [Peak],
    traces: &IonTraces,
    scan_dist: f32,
    ion_dist: f32,
    number_bins: u32,
) {
    struct Profile {
        peaks: Vec<usize>,
        last_row: usize,
        apex: u32,
        /// insertIonProfile puts the profile into bins `lo..hi`.
        lo: usize,
        hi: usize,
    }
    peaks.sort_by_key(|p| (p.row, p.first_scan));
    let nb = number_bins as usize;
    let dist = (number_bins as f32 * scan_dist).round() as usize;
    let len = nb * (traces.scan_count() + 1) + 1;
    // Instead of the old 2 * dist bin lists per profile, blocks of `dist` bins list every profile
    // whose range touches them (at most 3 blocks), in creation order like the old bin lists.
    let block = dist.max(1);
    let mut blocks: Vec<Vec<usize>> = vec![Vec::new(); len / block + 1];
    let mut profiles: Vec<Profile> = Vec::new();
    for k in 0..peaks.len() {
        let p = peaks[k];
        let b = (p.interpolated_max_scan() * nb as f64 + 0.5) as usize;
        if b >= len {
            continue; // the old code throws here; detected peaks never get this far
        }
        let next = |r: usize| (r + 1 < traces.len()).then_some(r + 1);
        let winner = blocks[b / block]
            .iter()
            .copied()
            .filter(|j| (profiles[*j].lo..profiles[*j].hi).contains(&b))
            .filter(|j| {
                let last = profiles[*j].last_row;
                next(last).is_some_and(|n| {
                    p.row as usize == n && (traces.mz[last] - traces.mz[n]).abs() <= ion_dist
                })
            })
            .min_by_key(|j| p.apex_scan.abs_diff(profiles[*j].apex));
        match winner {
            None => {
                // insertIonProfile: the new profile sits in bins [b - dist, b + dist).
                let (lo, hi) = (b.saturating_sub(dist), (b + dist).min(len - 1));
                if lo < hi {
                    (lo / block..=(hi - 1) / block).for_each(|i| blocks[i].push(profiles.len()));
                }
                profiles.push(Profile {
                    peaks: vec![k],
                    last_row: p.row as usize,
                    apex: summed_apex(traces, peaks, &[k], None).0,
                    lo,
                    hi,
                });
            }
            Some(j) if p.first_scan < profiles[j].apex && profiles[j].apex < p.last_scan => {
                let w = &mut profiles[j];
                w.peaks.push(k);
                w.last_row = p.row as usize;
                w.apex = summed_apex(traces, peaks, &w.peaks, None).0;
            }
            Some(_) => {}
        }
    }
    for w in profiles.iter().filter(|w| w.peaks.len() >= 2) {
        let (apex, precise) = summed_apex(traces, peaks, &w.peaks, Some(number_bins));
        for k in &w.peaks {
            let q = &mut peaks[*k];
            if (q.first_scan..=q.last_scan).contains(&apex) {
                q.apex_scan = apex;
                q.precise_apex = precise;
            }
        }
    }
}

/// Apex and precise apex of the old summed `Peak(list)` (§4.3): the members' intensities added
/// per scan over their joint range. The precise apex only with `bins` (the fit is the costly part).
fn summed_apex(
    traces: &IonTraces,
    peaks: &[Peak],
    members: &[usize],
    bins: Option<u32>,
) -> (u32, Option<f32>) {
    let first = members.iter().map(|k| peaks[*k].first_scan).min().unwrap();
    let last = members.iter().map(|k| peaks[*k].last_scan).max().unwrap();
    let mut sum = vec![0.0f32; (last - first + 1) as usize];
    for p in members.iter().map(|k| &peaks[*k]) {
        let trace = traces.trace(p.row as usize);
        for sc in p.first_scan..=p.last_scan {
            sum[(sc - first) as usize] += trace[sc as usize];
        }
    }
    let apex = argmax(&sum, 0);
    let precise = bins.and_then(|b| precise_apex(&sum, first as usize, apex, b));
    (first + apex as u32, precise)
}

/// Kernel of the old `tools::Filter::matched2ndDerivativeGaussianFilter` (§3.3): a negative second
/// derivative of a Gaussian, sigma = `density / 4`, cut at ±`density`.
fn matched_kernel(density: usize) -> Vec<f64> {
    let d = density.max(1) as isize;
    (-d..=d)
        .map(|x| {
            let u = 4.0 * x as f64 / d as f64;
            1.73 / (d as f64).sqrt() * (1.0 - u * u) * (-0.5 * u * u).exp()
        })
        .collect()
}

/// Correlation of `v` with `kernel` at index `i`; points outside the trace count as zero.
fn correlate<T: Copy + Into<f64>>(v: &[T], kernel: &[f64], i: usize) -> f64 {
    let d = kernel.len() / 2;
    let lo = i.saturating_sub(d);
    let hi = (i + d + 1).min(v.len());
    v[lo..hi]
        .iter()
        .zip(&kernel[lo + d - i..])
        .map(|(y, k)| k * (*y).into())
        .sum()
}

/// The old filter's last step: the first and last `d` values, where the kernel reaches past the
/// trace, take the value of the nearest scan that is not among them, or 0 if the trace is not longer
/// than `d` (§3.3).
fn hold_edges<T: Copy + Default>(f: &mut [T], d: usize) {
    let n = f.len();
    if d >= n {
        f.fill(T::default());
        return;
    }
    let first = f[d];
    f[..d].fill(first);
    let last = f[n - d - 1];
    f[n - d..].fill(last);
}

/// Old `tools::Filter::matched2ndDerivativeGaussianFilter` (§3.3): correlation with a negative
/// second derivative of a Gaussian, sigma = `density / 4`, cut at ±`density`, zeros outside the trace,
/// then the `density` values at each end held at the nearest inner value.
/// It responds most to Gaussian peaks with sigma ≈ `density / 4` scans.
pub fn matched_filter<T: Copy + Into<f64>>(v: &[T], density: usize) -> Vec<f64> {
    let kernel = matched_kernel(density);
    let mut f: Vec<f64> = (0..v.len()).map(|i| correlate(v, &kernel, i)).collect();
    hold_edges(&mut f, kernel.len() / 2);
    f
}

/// Old `MatchedPeakDetector<int,float>::detectPeaks` (§3.3–3.4) on one trace; the old
/// `MatchedPeakDetector(width, …)` uses density `4 * width`. Regions where the filtered trace rises
/// above `threshold`, more than 2 scans long. As in the old code, the peak metrics use the filtered
/// values with noise factor 0.5, there is no baseline check, and `signal_to_noise` keeps the old
/// default -1.
pub fn detect_peaks_matched(v: &[f32], density: usize, threshold: i32, bins: u32) -> Vec<Peak> {
    let (mut f, mut peaks) = (vec![0.0; v.len()], Vec::new());
    matched_peaks(
        v,
        &matched_kernel(density),
        threshold,
        bins,
        &mut f,
        &mut peaks,
    );
    peaks
}

/// Filters `v` into `f` (same length; truncated to integers like the old int filter) and appends
/// the matched-detector peaks of `f` to `out`.
fn matched_peaks(
    v: &[f32],
    kernel: &[f64],
    threshold: i32,
    bins: u32,
    f: &mut [f32],
    out: &mut Vec<Peak>,
) {
    for (i, y) in f.iter_mut().enumerate() {
        *y = correlate(v, kernel, i).trunc() as f32;
    }
    hold_edges(f, kernel.len() / 2);
    let thr = threshold as f32;
    let mut i = 1;
    while i < f.len() {
        // Starts where the previous scan lies below and this one above the threshold, ends at the
        // first scan at or below it (included); a region still open at the trace end is dropped.
        if f[i - 1] < thr && thr < f[i] {
            let Some(end) = (i + 1..f.len()).find(|j| f[*j] <= thr) else {
                break;
            };
            if end - i > 2 {
                let w = &f[i..=end];
                out.push(peak(w, i, argmax(w, 0), 0.5, -1.0, bins));
            }
            i = end;
        }
        i += 1;
    }
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

/// §3.6, old `determinePreciseMaximizationTime`: vertex of the least-squares parabola through the
/// values at apex − 1, apex, apex + 1 of `v = trace[first..]`, fitted in absolute bins (`bins` per
/// scan) like the old program and truncated to `1 / bins` scan after adding half a bin. `None` at
/// the peak edge or for a vertex ≥ 1 scan from the apex. The fit is the old one bit for bit
/// ([`crate::quadfit`]), so half-bin ties and flat tops land where the old program puts them.
fn precise_apex(v: &[f32], first: usize, apex: usize, bins: u32) -> Option<f32> {
    if apex == 0 || apex + 1 >= v.len() {
        return None;
    }
    let scan = (first + apex) as f64;
    let bins = bins as f64;
    let x = [-1.0, 0.0, 1.0].map(|k| bins * (scan + k));
    let y = [0, 1, 2].map(|k| v[apex - 1 + k] as f64);
    let c = crate::quadfit::quadratic(x, y);
    let q = c[1] / (c[2] * -2.0) + 0.5;
    // C's `(int)q`; outside the int range (or NaN) x86 gives INT_MIN, never within a scan.
    if q.is_nan() || q.abs() >= 2_147_483_648.0 {
        return None;
    }
    let t = q.trunc() / bins;
    ((t - scan).abs() < 1.0).then_some(t as f32)
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
    fn precise_apex_breaks_ties_like_the_old_fit() {
        // Vertex 0.25 scan right of the apex = a half-bin tie; the old fit rounds it up at scan 1
        // and 4 and down at scan 2 and 3 (its last bits depend on the absolute bins).
        let up = [1.0, 4.0, 3.0];
        assert_eq!(precise_apex(&up, 0, 1, 10), Some(1.3));
        assert_eq!(precise_apex(&up, 1, 1, 10), Some(2.2));
        assert_eq!(precise_apex(&up, 2, 1, 10), Some(3.2));
        assert_eq!(precise_apex(&up, 3, 1, 10), Some(4.3));
        assert_eq!(precise_apex(&[3.0, 4.0, 1.0], 0, 1, 10), Some(0.8));
        assert_eq!(precise_apex(&[3.0, 4.0, 1.0], 1, 1, 10), Some(1.7));
        assert_eq!(precise_apex(&up, 0, 1, 4), Some(1.25));
        // A flat top gives the rounding noise of the fit, as in the old program; zeros give NaN.
        assert_eq!(precise_apex(&[5.0, 5.0, 5.0], 0, 1, 10), Some(0.7));
        assert_eq!(precise_apex(&[0.0, 0.0, 0.0], 0, 1, 10), None);
        assert_eq!(precise_apex(&up, 0, 0, 10), None);
        assert_eq!(precise_apex(&up, 0, 2, 10), None);
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
        // Impulse response away from the 8 held scans at each end: the mirrored kernel.
        let mut impulse = [0.0; 33];
        impulse[16] = 1.0;
        let k = matched_filter(&impulse, 8);
        assert!(
            k[16] > 0.0 && k[8] < 0.0 && (k[8] - k[24]).abs() < 1e-12,
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
    fn matched_detection_filters_the_traces_in_place() {
        let mut v = noisy(300);
        gauss(&mut v, 150.0, 10.0, 50_000.0);
        let c = Chromatogram {
            scan_times_s: (0..v.len()).map(|s| s as f64).collect(),
            offsets: (0..=v.len()).collect(),
            mz: vec![73.0; v.len()],
            intensity: v.clone(),
            exact_masses: false,
            baseline: Default::default(),
        };
        let mut t = IonTraces::new(&c, 0.5).unwrap();
        let s = PeakSettings {
            detector: Detector::MATCHED,
            ..PeakSettings::default()
        };
        let p = t.detect_peaks(&s);
        assert_eq!(p, detect_peaks_matched(&v, 40, 3, 10), "{p:?}");
        assert!(p.iter().any(|p| p.apex_scan.abs_diff(150) <= 1), "{p:?}");
        let f: Vec<f32> = matched_filter(&v, 40)
            .iter()
            .map(|x| x.trunc() as f32)
            .collect();
        assert_eq!(t.trace(0), &f[..]);
        // The old filter holds the 40 scans at each end at the nearest inner value.
        let (head, tail) = (&f[..40], &f[260..]);
        assert!(head.iter().all(|y| *y == f[40]) && tail.iter().all(|y| *y == f[259]));
        assert_eq!(matched_filter(&[5.0f32; 3], 3), vec![0.0; 3]);
    }

    #[test]
    fn ion_traces_take_odd_mz_from_broken_files() {
        let c = Chromatogram {
            scan_times_s: vec![0.0],
            offsets: vec![0, 4],
            mz: vec![-0.5, 0.3, 0.5, f32::INFINITY],
            intensity: vec![1.0, 2.0, 4.0, 8.0],
            exact_masses: false,
            baseline: Default::default(),
        };
        let t = IonTraces::new(&c, 0.5).unwrap();
        assert_eq!(t.masses, vec![0, 1]);
        assert_eq!((t.trace(0), t.trace(1)), (&[3.0][..], &[4.0][..]));
    }

    #[test]
    fn centroiding_follows_the_old_cdf_import() {
        let mut out = vec![];
        // Close pairs: both points count when the smaller has half the maximum, else only the
        // larger (values seen in the old program's cache). 149.0 is its own peak.
        centroid_scan(
            &[147.047, 147.089, 148.047, 148.089, 149.0],
            &[16012.0, 15171.0, 140.0, 339.0, 7.0],
            0.5,
            &mut out,
        );
        assert_eq!(out.len(), 3);
        assert_eq!((out[0].1, out[1].1, out[2]), (31183.0, 339.0, (149.0, 7.0)));
        assert!((out[0].0 - 147.0674).abs() < 1e-3 && (out[1].0 - 148.089).abs() < 1e-4);
        // A valley (two falls, then two rises) splits points closer than 0.1.
        out.clear();
        let m: Vec<f32> = (0..7).map(|k| 100.0 + 0.01 * k as f32).collect();
        centroid_scan(&m, &[5.0, 3.0, 1.0, 2.0, 4.0, 6.0, 1.0], 0.5, &mut out);
        assert_eq!(out.iter().map(|p| p.1).collect::<Vec<_>>(), vec![8.0, 10.0]);
        // A higher SF_FWHM drops the smaller point of the pair: mass and intensity of the top alone.
        out.clear();
        centroid_scan(&[147.047, 147.089], &[16012.0, 15171.0], 0.99, &mut out);
        assert_eq!(out.len(), 1);
        assert!(
            (out[0].0 - 147.047).abs() < 1e-4 && out[0].1 == 16012.0,
            "{out:?}"
        );
        // Close points of one nominal mass merge in the nominal traces; others stay as they are.
        let c = Chromatogram {
            scan_times_s: vec![0.0],
            offsets: vec![0, 3],
            mz: vec![73.047, 73.089, 74.0],
            intensity: vec![10.0, 2.0, 5.0],
            exact_masses: false,
            baseline: Default::default(),
        };
        let t = IonTraces::new(&c, 0.5).unwrap();
        assert_eq!((t.trace(0), t.trace(1)), (&[10.0][..], &[5.0][..]));
    }

    #[test]
    fn reorder_peaks_aligns_apexes_of_split_masses() {
        // 73.047 and 73.089 (0.042 apart) carry one fragment with apexes 3 scans apart; 74.0 is a
        // whole mass away and keeps its apex.
        let scans = 200;
        let mut traces = [noisy(scans), noisy(scans), noisy(scans)];
        gauss(&mut traces[0], 100.0, 3.0, 20_000.0);
        gauss(&mut traces[1], 103.0, 3.0, 10_000.0);
        gauss(&mut traces[2], 103.0, 3.0, 10_000.0);
        let c = Chromatogram {
            scan_times_s: (0..scans).map(|s| s as f64).collect(),
            offsets: (0..=scans).map(|s| 3 * s).collect(),
            mz: (0..scans).flat_map(|_| [73.047f32, 73.089, 74.0]).collect(),
            intensity: (0..scans)
                .flat_map(|s| [traces[0][s], traces[1][s], traces[2][s]])
                .collect(),
            exact_masses: false,
            baseline: Default::default(),
        };
        let mut t = IonTraces::exact(&c).unwrap();
        let mut p = t.detect_peaks(&PeakSettings::default());
        let apexes = |p: &[Peak]| p.iter().map(|p| p.apex_scan).collect::<Vec<_>>();
        assert_eq!(apexes(&p), vec![100, 103, 103]);
        reorder_peaks(&mut p, &t, 10.0, 0.05, 10);
        // The summed profile of rows 0 and 1 peaks at scan 101.
        assert_eq!(apexes(&p), vec![101, 101, 103]);
        assert_eq!(p[0].precise_apex, p[1].precise_apex);
        // Apexes further apart than scan_dist stay apart.
        let mut q = t.detect_peaks(&PeakSettings::default());
        reorder_peaks(&mut q, &t, 0.2, 0.05, 10);
        assert_eq!(apexes(&q), vec![100, 103, 103]);
        // detectAllPeaks(true, 10, 0.05) reorders exact traces only.
        assert_eq!(t.detect_all_peaks(&PeakSettings::default()), p);
        let mut nominal = IonTraces::new(&c, 0.5).unwrap();
        assert_eq!(
            nominal.detect_all_peaks(&PeakSettings::default()),
            nominal.detect_peaks(&PeakSettings::default())
        );
    }

    #[test]
    fn reorder_peaks_gives_a_tie_to_the_older_profile() {
        // Row 0 has two peaks (95 and 105); the row-1 peak at 100 is 5 scans from both. The old
        // bin lists are appended and walked from the front, so the profile made first (95) wins.
        let scans = 300;
        let mut traces = [noisy(scans), noisy(scans)];
        gauss(&mut traces[0], 95.0, 2.0, 20_000.0);
        gauss(&mut traces[0], 105.0, 2.0, 20_000.0);
        gauss(&mut traces[1], 100.0, 3.0, 20_000.0);
        let c = Chromatogram {
            scan_times_s: (0..scans).map(|s| s as f64).collect(),
            offsets: (0..=scans).map(|s| 2 * s).collect(),
            mz: (0..scans).flat_map(|_| [73.047f32, 73.089]).collect(),
            intensity: (0..scans)
                .flat_map(|s| [traces[0][s], traces[1][s]])
                .collect(),
            exact_masses: false,
            baseline: Default::default(),
        };
        let t = IonTraces::exact(&c).unwrap();
        let peak = |row, first, apex, last| Peak {
            row,
            first_scan: first,
            apex_scan: apex,
            last_scan: last,
            precise_apex: Some(apex as f32),
            signal_to_noise: 10.0,
            quality: 0.0,
            sharpness: 1.0,
        };
        let mut p = [
            peak(0, 80, 95, 100),
            peak(0, 101, 105, 120),
            peak(1, 80, 100, 120),
        ];
        reorder_peaks(&mut p, &t, 10.0, 0.05, 10);
        assert_eq!(p[0].apex_scan, p[2].apex_scan);
        assert!(p[2].apex_scan < 100);
        assert_eq!(p[1].apex_scan, 105);
    }

    #[test]
    fn reorder_peaks_leaves_out_the_last_bin_of_a_profile() {
        // The profile of the first peak (bin 1003) covers bins 903..1103 (insertIonProfile, §4.2); bin
        // 1103 shares a lookup block with the covered bins, so the range test decides.
        let scans = 300;
        let mut traces = [noisy(scans), noisy(scans)];
        gauss(&mut traces[0], 100.0, 3.0, 20_000.0);
        gauss(&mut traces[1], 105.0, 3.0, 20_000.0);
        let c = Chromatogram {
            scan_times_s: (0..scans).map(|s| s as f64).collect(),
            offsets: (0..=scans).map(|s| 2 * s).collect(),
            mz: (0..scans).flat_map(|_| [73.047f32, 73.089]).collect(),
            intensity: (0..scans)
                .flat_map(|s| [traces[0][s], traces[1][s]])
                .collect(),
            exact_masses: false,
            baseline: Default::default(),
        };
        let t = IonTraces::exact(&c).unwrap();
        let peak = |row, apex, precise| Peak {
            row,
            first_scan: 90,
            apex_scan: apex,
            last_scan: 120,
            precise_apex: Some(precise),
            signal_to_noise: 10.0,
            quality: 0.0,
            sharpness: 1.0,
        };
        let mut joined = [peak(0, 100, 100.3), peak(1, 105, 110.2)];
        reorder_peaks(&mut joined, &t, 10.0, 0.05, 10);
        assert_eq!(joined[0].apex_scan, joined[1].apex_scan);
        let mut apart = [peak(0, 100, 100.3), peak(1, 105, 110.3)];
        reorder_peaks(&mut apart, &t, 10.0, 0.05, 10);
        assert_eq!((apart[0].apex_scan, apart[1].apex_scan), (100, 105));
    }

    #[test]
    fn exact_traces_keep_close_masses_apart() {
        // Scan 1 lacks 73.089, which reads 0 there.
        let c = Chromatogram {
            scan_times_s: vec![0.0, 1.0],
            offsets: vec![0, 3, 5],
            mz: vec![73.047, 73.089, 147.0, 73.047, 147.0],
            intensity: vec![1.0, 2.0, 3.0, 4.0, 5.0],
            exact_masses: false,
            baseline: Default::default(),
        };
        assert_eq!(IonTraces::new(&c, 0.5).unwrap().trace(0), &[3.0, 4.0][..]);
        let t = IonTraces::exact(&c).unwrap();
        assert_eq!(t.masses, vec![73, 73, 147]);
        assert_eq!(t.mz, vec![73.047, 73.089, 147.0]);
        assert_eq!(t.trace(0), &[1.0, 4.0][..]);
        assert_eq!(t.trace(1), &[2.0, 0.0][..]);
        assert_eq!(t.trace(2), &[3.0, 5.0][..]);
        // A FastFlight run gets the exact traces from `new` too.
        let ff = Chromatogram {
            exact_masses: true,
            ..c
        };
        assert_eq!(IonTraces::new(&ff, 0.5).unwrap().mz, t.mz);
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
            exact_masses: false,
            baseline: Default::default(),
        };
        let mut t = IonTraces::new(&c, 0.5).unwrap();
        assert_eq!(t.masses, vec![73, 147]);
        assert_eq!(t.mz, vec![73.0, 147.0]);
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
            exact_masses: false,
            baseline: Default::default(),
        };
        assert!(
            IonTraces::new(&empty, 0.5)
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
            exact_masses: false,
            baseline: Default::default(),
        };
        assert_eq!(
            IonTraces::new(&c, 0.5),
            Err(TooLarge {
                masses: 65_536,
                scans
            })
        );
    }
}
