//! Deconvolution: groups the ion-trace peaks into components with one spectrum each, the old
//! `GCMSAnalyzer::deconvoluteChromatogram` (docs/reverse-engineering/algorithms.md §4.4–4.9).
//! Deviations are listed in docs/architecture.md.

use std::collections::BTreeMap;

use rayon::prelude::*;

use crate::peaks::{IonTraces, Peak, matched_filter};

/// The old `AN_*` deconvolution settings; defaults are the old static values (file-formats.md §5.1).
#[derive(Debug, Clone, PartialEq)]
pub struct DeconvSettings {
    /// `AN_DECONVOLUTION_WIDTH`: width in scans of the matched filter over the apex histogram.
    pub width: f64,
    /// `AN_NUMBER_BINS`: histogram bins per scan; keep it equal to `PeakSettings::number_bins`.
    pub number_bins: u32,
    /// `AN_MIN_PEAK_NUMBER`: minimum ion peaks per component, before and after the fit.
    pub min_peak_number: usize,
    /// `AN_MIN_HEIGHT_BP`: ions below this fraction of the component's base peak leave its spectrum.
    pub min_height_bp: f64,
    /// `AN_MODELED_VALUES_FOR_INTEGRATION`: integrate the fitted model instead of the raw trace.
    pub modeled_integration: bool,
    /// `QUANT_EXCL_IONS`: inclusive mass ranges never used to build the model peak.
    pub excluded_model_ions: Vec<(f32, f32)>,
}

impl Default for DeconvSettings {
    fn default() -> Self {
        Self {
            width: 1.0,
            number_bins: 10,
            min_peak_number: 10,
            min_height_bp: 0.0,
            modeled_integration: true,
            excluded_model_ions: Vec::new(),
        }
    }
}

/// A deconvoluted substance: the old `Compound` before RI and library search.
#[derive(Debug, Clone, PartialEq)]
pub struct Component {
    /// Scan range and first highest scan of the model peak.
    pub first_scan: u32,
    pub apex_scan: u32,
    pub last_scan: u32,
    /// Model peak: summed raw intensities of the model peaks for each scan in `first_scan..=last_scan`.
    pub model: Vec<f64>,
    /// Deconvoluted spectrum as (nominal mass, modelled apex height), ascending mass.
    pub spectrum: Vec<(u32, u32)>,
    /// Ion peaks behind the spectrum, sharpest first. Like the old code a peak can appear twice (§4.5).
    pub peaks: Vec<FittedPeak>,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FittedPeak {
    /// Index into the peak list given to `deconvolute`.
    pub peak: usize,
    /// Share of the model peak in this ion trace (`c2` of the fit).
    pub model_share: f64,
    /// Plain sum over the overlap with the model of the fitted (or raw) intensities.
    pub integral: f64,
}

/// Old `deconvoluteChromatogram` (§4.6) for the peaks of `traces` (from `IonTraces::detect_peaks`).
/// Builds the apex histogram once, then fits the components in parallel.
pub fn deconvolute(traces: &IonTraces, peaks: &[Peak], s: &DeconvSettings) -> Vec<Component> {
    let nb = s.number_bins.max(1) as usize;
    let len = nb * (traces.scan_count() + 1) + 1;
    // §4.4: apex intensities summed per sub-scan bin, with the peaks of each bin (counting sort keeps
    // the old order: by trace, then scan). Without a precise apex the old field keeps the apex scan.
    let bin_of: Vec<Option<usize>> = peaks
        .iter()
        .map(|p| {
            let t = p.precise_apex.unwrap_or(p.apex_scan as f32) as f64;
            let b = (nb as f64 * t + 0.5) as usize;
            (b < len).then_some(b)
        })
        .collect();
    let mut start = vec![0usize; len + 1];
    bin_of.iter().flatten().for_each(|b| start[b + 1] += 1);
    (0..len).for_each(|b| start[b + 1] += start[b]);
    let (mut fill, mut members, mut sum) = (start.clone(), vec![0; start[len]], vec![0.0f64; len]);
    for (k, b) in bin_of.iter().enumerate() {
        if let Some(b) = *b {
            members[fill[b]] = k;
            fill[b] += 1;
            let p = &peaks[k];
            sum[b] += traces.trace(p.row as usize)[p.apex_scan as usize] as f64;
        }
    }
    let f = matched_filter(&sum, (nb as f64 * s.width) as usize);
    regions(&f, &start, &members)
        .into_par_iter()
        .filter_map(|t| fit_component(traces, peaks, t, s))
        .collect()
}

/// §4.5: regions where the filtered apex histogram (truncated to integers) rises above 5, plus the
/// falling flanks on both sides. Keeps the old quirks: the flank walk adds the first and last bin
/// of the region a second time, and a region of one bin is dropped.
fn regions(f: &[f64], start: &[usize], members: &[usize]) -> Vec<Vec<usize>> {
    const T: f64 = 5.0;
    let add = |t: &mut Vec<usize>, b: usize| t.extend_from_slice(&members[start[b]..start[b + 1]]);
    let mut out = Vec::new();
    let Some(f0) = f.first() else { return out };
    let mut prev = f0.trunc();
    let mut i = 0;
    while i < f.len() {
        let (begin, mut end) = (i, i);
        let mut t = Vec::new();
        if prev < T && f[i].trunc() > T {
            add(&mut t, i);
            loop {
                i += 1;
                if i >= f.len() {
                    break;
                }
                add(&mut t, i);
                end = i;
                if f[i].trunc() <= T {
                    break;
                }
            }
        }
        if begin != end {
            let mut n = begin;
            while n > 1 && f[n] > f[n - 1] {
                add(&mut t, n);
                n -= 1;
            }
            let mut n = end;
            while n + 1 < f.len() && f[n] > f[n + 1] {
                add(&mut t, n);
                n += 1;
            }
            out.push(t);
        }
        if let Some(x) = f.get(i) {
            prev = x.trunc();
        }
        i += 1;
    }
    out
}

/// §4.6–4.9 for one region: build the model peak from the sharpest clean peaks, fit every ion
/// trace against it and keep the modelled apex heights as the spectrum.
fn fit_component(
    traces: &IonTraces,
    peaks: &[Peak],
    mut t: Vec<usize>,
    s: &DeconvSettings,
) -> Option<Component> {
    if t.is_empty() || t.len() < s.min_peak_number {
        return None;
    }
    t.sort_by(|a, b| peaks[*b].sharpness.total_cmp(&peaks[*a].sharpness));
    let mass = |k: usize| traces.masses[peaks[k].row as usize];
    let excluded = |m: u32| {
        s.excluded_model_ions
            .iter()
            .any(|(lo, hi)| *lo <= m as f32 && m as f32 <= *hi)
    };
    let quota = ((t.len() as f64 * 0.25) as usize).max(1);
    let mut m: Vec<usize> = t
        .iter()
        .copied()
        .filter(|k| peaks[*k].quality <= 10.0 && !excluded(mass(*k)))
        .take(quota)
        .collect();
    if m.is_empty() {
        let s0 = peaks[t[0]].sharpness as f64 * 0.75;
        m = t
            .iter()
            .copied()
            .take_while(|k| s0 <= peaks[*k].sharpness as f64)
            .collect();
    }

    let first = m.iter().map(|k| peaks[*k].first_scan).min()?;
    let last = m.iter().map(|k| peaks[*k].last_scan).max()?;
    let mut model = vec![0.0f64; (last - first + 1) as usize];
    for k in &m {
        let p = &peaks[*k];
        let trace = traces.trace(p.row as usize);
        for sc in p.first_scan..=p.last_scan {
            model[(sc - first) as usize] += trace[sc as usize] as f64;
        }
    }
    let apex = (0..model.len()).fold(0, |a, k| if model[k] > model[a] { k } else { a });
    let model_max = model[apex];
    let at = |sc: u32| model[(sc - first) as usize];

    let mut fitted = Vec::new();
    for &k in &t {
        let p = &peaks[k];
        let (lo, hi) = (first.max(p.first_scan), last.min(p.last_scan));
        if hi < lo + 2 {
            continue;
        }
        let trace = traces.trace(p.row as usize);
        let y = |sc: u32| trace[sc as usize] as f64;
        let Some(c2) = model_share(lo..=hi, at, y) else {
            continue;
        };
        let c2 = round_sig(c2, 8);
        let h = (c2 * model_max).trunc();
        if h < 1.0 {
            continue;
        }
        let integral = (lo..=hi)
            .map(|sc| {
                if s.modeled_integration {
                    c2 * at(sc)
                } else {
                    y(sc)
                }
            })
            .sum();
        fitted.push((k, c2, h as u32, integral));
    }
    if fitted.len() < s.min_peak_number {
        return None;
    }

    // A mass seen twice keeps its first (sharpest) height; the base peak counts all of them.
    let bp = fitted.iter().map(|f| f.2).max().unwrap_or(0);
    let mut spectrum = BTreeMap::new();
    for f in &fitted {
        spectrum.entry(mass(f.0)).or_insert(f.2);
    }
    spectrum.retain(|_, h| *h as f64 >= s.min_height_bp * bp as f64);
    Some(Component {
        first_scan: first,
        apex_scan: first + apex as u32,
        last_scan: last,
        peaks: fitted
            .iter()
            .filter(|f| spectrum.contains_key(&mass(f.0)))
            .map(|&(peak, model_share, _, integral)| FittedPeak {
                peak,
                model_share,
                integral,
            })
            .collect(),
        spectrum: spectrum.into_iter().collect(),
        model,
    })
}

/// §4.8: `c2` of the least-squares fit `y(scan) ≈ c0 + c1 * scan + c2 * model(scan)`.
/// Centred normal equations; `None` where model and scan are collinear.
// shortcut: the old gsl SVD returns a minimum-norm fit for collinear cases, here the peak is dropped;
// port the SVD if the golden test shows such peaks.
fn model_share(
    scans: std::ops::RangeInclusive<u32>,
    model: impl Fn(u32) -> f64,
    y: impl Fn(u32) -> f64,
) -> Option<f64> {
    let n = scans.clone().count() as f64;
    let (sx, sm, sy) = scans.clone().fold((0.0, 0.0, 0.0), |(a, b, c), sc| {
        (a + sc as f64, b + model(sc), c + y(sc))
    });
    let (mx, mm, my) = (sx / n, sm / n, sy / n);
    let [mut xx, mut xm, mut mmm, mut xy, mut my2] = [0.0f64; 5];
    for sc in scans {
        let (dx, dm, dy) = (sc as f64 - mx, model(sc) - mm, y(sc) - my);
        xx += dx * dx;
        xm += dx * dm;
        mmm += dm * dm;
        xy += dx * dy;
        my2 += dm * dy;
    }
    let det = xx * mmm - xm * xm;
    (det > 1e-12 * xx * mmm).then(|| (xx * my2 - xm * xy) / det)
}

/// The old `PeakEssential::roundSignificantFigures`.
fn round_sig(x: f64, digits: i32) -> f64 {
    if x == 0.0 || !x.is_finite() {
        return x;
    }
    let scale = 10f64.powi(digits - 1 - x.abs().log10().floor() as i32);
    (x * scale).round() / scale
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Chromatogram;
    use crate::peaks::PeakSettings;

    #[test]
    fn regions_keep_the_old_flank_quirks() {
        // One peak per bin, numbered like the bin.
        let f = [0.0, 0.0, 3.0, 10.0, 20.0, 10.0, 3.0, 0.0, 0.0];
        let start: Vec<usize> = (0..=f.len()).collect();
        let members: Vec<usize> = (0..f.len()).collect();
        assert_eq!(
            regions(&f, &start, &members),
            vec![vec![3, 4, 5, 6, 3, 2, 6]]
        );
        assert!(regions(&[0.0, 9.0], &start, &members).is_empty());
        assert!(regions(&[], &start, &members).is_empty());
    }

    #[test]
    fn model_share_recovers_scaled_model_over_a_baseline() {
        let model = |sc: u32| (-0.5 * ((sc as f64 - 50.0) / 2.0).powi(2)).exp() * 1000.0;
        let c2 = model_share(40..=60, model, |sc| {
            0.3 * model(sc) + 20.0 + 0.5 * sc as f64
        });
        assert!((c2.unwrap() - 0.3).abs() < 1e-9);
        assert_eq!(model_share(40..=60, |_| 1.0, |_| 1.0), None);
        assert_eq!(round_sig(0.123456789, 8), 0.12345679);
    }

    /// Apex scan and (mass, height) of each ion.
    type Substance = (f64, Vec<(u32, f64)>);

    /// Two substances 6 scans apart, each with its own 12 masses, over a noisy baseline.
    fn two_substances() -> (Chromatogram, [Substance; 2]) {
        let subst = |apex: f64, masses: std::ops::Range<u32>| {
            let ions: Vec<(u32, f64)> = masses
                .map(|m| (m, 3_000.0 + 4_000.0 * (m % 7) as f64))
                .collect();
            (apex, ions)
        };
        let s = [subst(100.0, 50..62), subst(106.0, 62..74)];
        let scans = 300u32;
        let mut x = 12345u32;
        let mut noise = move || {
            x = x.wrapping_mul(1_103_515_245).wrapping_add(12345);
            ((x >> 16) % 41) as f32 - 20.0
        };
        let mut intensity = Vec::new();
        for sc in 0..scans {
            for m in 50..74u32 {
                let (apex, ions) = &s[usize::from(m >= 62)];
                let h = ions.iter().find(|(k, _)| *k == m).unwrap().1;
                let g = h * (-0.5 * ((sc as f64 - apex) / 2.5).powi(2)).exp();
                intensity.push(1000.0 + noise() + g as f32);
            }
        }
        let c = Chromatogram {
            scan_times_s: (0..scans).map(f64::from).collect(),
            offsets: (0..=scans as usize).map(|s| 24 * s).collect(),
            mz: (0..scans)
                .flat_map(|_| (50..74).map(|m| m as f32))
                .collect(),
            intensity,
        };
        (c, s)
    }

    #[test]
    fn separates_two_overlapping_substances() {
        let (c, truth) = two_substances();
        let traces = IonTraces::new(&c);
        let peaks = traces.detect_peaks(&PeakSettings::default());
        let comps = deconvolute(&traces, &peaks, &DeconvSettings::default());
        assert_eq!(comps.len(), 2, "{comps:?}");
        for (comp, (apex, ions)) in comps.iter().zip(&truth) {
            assert_eq!(comp.apex_scan as f64, *apex);
            let masses: Vec<u32> = comp.spectrum.iter().map(|(m, _)| *m).collect();
            assert_eq!(masses, ions.iter().map(|(m, _)| *m).collect::<Vec<_>>());
            // The model sums raw intensities, baseline included, so heights share one scale factor;
            // relative to the base peak they match the true ones.
            let bp = comp.spectrum.iter().map(|(_, h)| *h).max().unwrap() as f64;
            let want_bp = ions.iter().map(|(_, h)| *h).fold(0.0, f64::max);
            for ((_, h), (_, want)) in comp.spectrum.iter().zip(ions) {
                let (got, want) = (*h as f64 / bp, want / want_bp);
                assert!((got - want).abs() < 0.03, "{got} vs {want}");
            }
            assert!(comp.peaks.iter().all(|p| p.integral > 0.0));
        }
    }

    #[test]
    fn rejects_components_with_too_few_ions_and_empty_input() {
        let (c, _) = two_substances();
        let traces = IonTraces::new(&c);
        let peaks = traces.detect_peaks(&PeakSettings::default());
        let strict = DeconvSettings {
            min_peak_number: 13,
            ..Default::default()
        };
        assert!(deconvolute(&traces, &peaks, &strict).is_empty());
        assert!(deconvolute(&traces, &[], &DeconvSettings::default()).is_empty());
    }
}
