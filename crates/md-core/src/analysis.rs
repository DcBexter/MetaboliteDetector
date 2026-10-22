//! Analyze run: ion traces, peaks and deconvolution of one run behind one interface, for the app,
//! the batch and its rescan.

use std::borrow::Borrow;

use crate::Chromatogram;
use crate::deconv::{Component, DeconvSettings, deconvolute};
use crate::job::Cancel;
use crate::peaks::{IonTraces, Peak, PeakSettings, TooLarge};

/// Which peak detection a caller wants (algorithms.md §4.0).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PeakPass {
    /// Main window and batch: old `detectAllPeaks(reorder=true)`, [`IonTraces::detect_all_peaks`].
    Reordered,
    /// Batch rescan: old `ResultSeries::scanIons` with a fresh `DefaultPeakDetector`,
    /// [`IonTraces::detect_peaks`] without `reorder_peaks`.
    Fresh,
}

/// The peak stage of a run: the ion-major index, the scan times and the peaks.
pub struct RunPeaks {
    pub traces: IonTraces,
    pub times: Vec<f64>,
    pub peaks: Vec<Peak>,
    /// [`PeakPass::Fresh`]'s peaks (`peaks` before `reorder_peaks`), whatever the pass: the batch
    /// rescan reuses them from pass 1.
    pub fresh: Vec<Peak>,
}

/// A deconvoluted run: components in detection order, their peaks and each peak's trace mass.
pub struct RunAnalysis {
    pub times: Vec<f64>,
    pub peaks: Vec<Peak>,
    pub components: Vec<Component>,
    /// Mass of the trace of each peak (`peaks[i].row`).
    pub masses: Vec<u32>,
    /// [`RunPeaks::fresh`] and the trace mass of each.
    pub fresh: Vec<Peak>,
    pub fresh_masses: Vec<u32>,
}

/// Why a run was not analysed.
#[derive(Debug, Clone, Copy, thiserror::Error, PartialEq, Eq)]
pub enum AnalyzeError {
    #[error(transparent)]
    TooLarge(#[from] TooLarge),
    #[error("cancelled")]
    Cancelled,
}

/// Builds the ion traces and detects the peaks. An owned `run` is dropped right after the index is
/// built: the spectra are not needed any more, which halves the peak memory per run.
pub fn find_run_peaks(
    run: impl Borrow<Chromatogram>,
    s: &PeakSettings,
    pass: PeakPass,
    cancel: &Cancel,
) -> Result<RunPeaks, AnalyzeError> {
    let mut traces = IonTraces::new(run.borrow(), s.fwhm)?;
    let times = run.borrow().scan_times_s.clone();
    drop(run);
    // The stages are not interruptible inside; `cancel` is looked at between them.
    if cancel.is_cancelled() {
        return Err(AnalyzeError::Cancelled);
    }
    let fresh = traces.detect_peaks(s);
    let mut peaks = fresh.clone();
    if pass == PeakPass::Reordered {
        traces.reorder(&mut peaks, s);
    }
    Ok(RunPeaks {
        traces,
        times,
        peaks,
        fresh,
    })
}

/// [`find_run_peaks`], then the deconvolution of all peaks; the ion index is dropped. A set `cancel`
/// gives [`AnalyzeError::Cancelled`] at the next stage boundary.
pub fn analyze_run(
    run: impl Borrow<Chromatogram>,
    peak_settings: &PeakSettings,
    deconv: &DeconvSettings,
    pass: PeakPass,
    cancel: &Cancel,
) -> Result<RunAnalysis, AnalyzeError> {
    let RunPeaks {
        traces,
        times,
        peaks,
        fresh,
    } = find_run_peaks(run, peak_settings, pass, cancel)?;
    if cancel.is_cancelled() {
        return Err(AnalyzeError::Cancelled);
    }
    let components = deconvolute(&traces, &peaks, &times, deconv);
    let mass_of = |p: &[Peak]| p.iter().map(|p| traces.masses[p.row as usize]).collect();
    Ok(RunAnalysis {
        masses: mass_of(&peaks),
        fresh_masses: mass_of(&fresh),
        times,
        peaks,
        components,
        fresh,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    /// An exact-mass run: 73.047 and 73.089 carry one fragment with apexes 3 scans apart (merged by
    /// `reorder_peaks`), 74.0 is a whole mass away.
    fn split_masses() -> Chromatogram {
        let scans = 200;
        let mut x = 12345u32;
        let mut trace = |center: f64, height: f64| -> Vec<f32> {
            (0..scans)
                .map(|k| {
                    x = x.wrapping_mul(1_103_515_245).wrapping_add(12345);
                    let noise = 1000.0 + ((x >> 16) % 41) as f64 - 20.0;
                    (noise + height * (-0.5 * ((k as f64 - center) / 3.0).powi(2)).exp()) as f32
                })
                .collect()
        };
        let t = [
            trace(100.0, 20_000.0),
            trace(103.0, 10_000.0),
            trace(103.0, 10_000.0),
        ];
        Chromatogram {
            scan_times_s: (0..scans).map(|s| s as f64).collect(),
            offsets: (0..=scans).map(|s| 3 * s).collect(),
            mz: (0..scans).flat_map(|_| [73.047f32, 73.089, 74.0]).collect(),
            intensity: (0..scans)
                .flat_map(|s| [t[0][s], t[1][s], t[2][s]])
                .collect(),
            exact_masses: true,
            baseline: Default::default(),
        }
    }

    #[test]
    fn peak_pass_selects_reordering() {
        let (c, s) = (split_masses(), PeakSettings::default());
        let apexes = |pass| -> Vec<u32> {
            let p = find_run_peaks(&c, &s, pass, &Cancel::default())
                .unwrap()
                .peaks;
            p.iter().map(|p| p.apex_scan).collect()
        };
        assert_eq!(apexes(PeakPass::Fresh), vec![100, 103, 103]);
        assert_eq!(apexes(PeakPass::Reordered), vec![101, 101, 103]);
    }

    #[test]
    /// Pins the refactor: `analyze_run` gives what the callers built by hand before #361.
    fn analyze_run_is_the_chain_for_both_passes() {
        let (c, s, d) = (
            split_masses(),
            PeakSettings::default(),
            DeconvSettings {
                min_peak_number: 1,
                ..Default::default()
            },
        );
        let fresh = (IonTraces::new(&c, s.fwhm).unwrap()).detect_peaks(&s);
        for pass in [PeakPass::Reordered, PeakPass::Fresh] {
            let mut traces = IonTraces::new(&c, s.fwhm).unwrap();
            let peaks = match pass {
                PeakPass::Reordered => traces.detect_all_peaks(&s),
                PeakPass::Fresh => traces.detect_peaks(&s),
            };
            let components = deconvolute(&traces, &peaks, &c.scan_times_s, &d);
            assert!(!components.is_empty());
            let masses: Vec<u32> = (peaks.iter())
                .map(|p| traces.masses[p.row as usize])
                .collect();
            // Borrowed (app) and owned (batch) runs give the same result.
            for a in [
                analyze_run(&c, &s, &d, pass, &Cancel::default()).unwrap(),
                analyze_run(c.clone(), &s, &d, pass, &Cancel::default()).unwrap(),
            ] {
                assert_eq!(a.peaks, peaks);
                assert_eq!(a.components, components);
                assert_eq!(a.times, c.scan_times_s);
                assert_eq!(a.masses, masses);
                // The rescan's list, whatever the pass.
                assert_eq!(a.fresh, fresh);
                assert_eq!(a.fresh_masses.len(), fresh.len());
            }
        }
        // The pass reaches the result: the reordered apexes differ from the fresh ones.
        let analyzed = |pass| {
            analyze_run(&c, &s, &d, pass, &Cancel::default())
                .unwrap()
                .peaks
        };
        assert_ne!(analyzed(PeakPass::Reordered), analyzed(PeakPass::Fresh));
    }

    #[test]
    fn a_set_cancel_stops_before_the_peaks() {
        let (c, s, d) = (
            split_masses(),
            PeakSettings::default(),
            DeconvSettings::default(),
        );
        let cancel = Cancel::default();
        cancel.cancel();
        let p = PeakPass::Reordered;
        assert_eq!(
            find_run_peaks(&c, &s, p, &cancel).err(),
            Some(AnalyzeError::Cancelled)
        );
        assert_eq!(
            analyze_run(&c, &s, &d, p, &cancel).err(),
            Some(AnalyzeError::Cancelled)
        );
    }
}
