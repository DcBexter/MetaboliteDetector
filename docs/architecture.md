# Architecture

MetaboliteDetector is being rebuilt in Rust as a native desktop app for Windows, Linux and macOS.
The old Qt5/C++ version (3.3, 2020) is the behavioural reference; its features are listed in
[reverse-engineering/feature-inventory.md](reverse-engineering/feature-inventory.md).

## Decisions

| Topic | Decision | Why |
|-------|----------|-----|
| Language | Rust (edition 2024) | Memory-safe, fast numerics, one static binary per OS, no runtime to install. |
| GUI | [egui](https://github.com/emilk/egui)/eframe + egui_plot | Pure Rust, single executable, interactive plots that handle large series. Chosen by the maintainer 2026-10-19. |
| File dialogs | rfd | Native dialogs on every OS. |
| UI language | `sys-locale`, app only | The system language as the default (S8); already in the tree through winit. See [UI language](#ui-language). |
| NetCDF | `netcdf3` (pure Rust) | ANDI-MS files are NetCDF-3; no libnetcdf/HDF5 C dependency to ship on Windows. |
| Parallelism | rayon | Data-parallel loops over scans, ions and files. |
| Updates | `ureq` (rustls), `serde_json`, `sha2`, `semver`, `tar` + `flate2`, app only | `ureq`: HTTPS to the GitHub API and release assets without a system TLS library; `serde_json`: read the release list; `sha2`: verify downloads against `SHA256SUMS.txt`; `semver`: order versions incl. pre-releases; `tar`: unpack the Linux tar.gz. See [Updates](#updates). |
| R console | Not ported directly | RInside tied the old app to a Linux R install. Decision (#26): no embedded R. The batch table exports as CSV (`Table::write_csv`, shared by app and md-cli) and `scripts/r/load_results.R` loads it as `data` like the old console. Optional embedded scripting later. |
| Search cutoff | LS_CUTOFF everywhere, also the main window | one rule; hits below are not kept, lowering the cutoff re-runs the identification; applies to the NIST list too (#399) |

## Layout

```
Cargo.toml              workspace: members crates/*, shared dependencies and lints, release profile
rust-toolchain.toml     pinned Rust version with clippy and rustfmt; Renovate raises it
crates/md-core/         library: data model, import/export, analysis algorithms (no GUI deps)
  src/lib.rs            Chromatogram (flat ANDI-style arrays), TIC/SIC, plot decimation
  src/andi.rs           ANDI-MS NetCDF import with input validation, NetCDF export
  src/diskscan.rs       old raw-data cache (.bin + .idx) reader, RI calibration stored in the .idx
  src/library.rs        compound libraries: old .lbr (versions 2-5) reader, AMDIS MSL import/export,
                        library text format (.mdlib) read/write
  src/peaks.rs          ion traces (ion-major index), noise, baseline and peak detection per trace
  src/deconv.rs         deconvolution of ion-trace peaks into components with spectra
  src/analysis.rs       analyze run: ion traces, peaks (PeakPass) and deconvolution of one run, the one
                        chain behind the app, the batch and its rescan
  src/settings.rs       old GCMSSettings: values and defaults, `analysis()` (the one derivation of AnalysisSettings),
                        `KEYS` with the key names in `settings::keys`, old settings file read and write
  src/job.rs            `Cancel` (cooperative stop flag) and `Progress` (done/total steps) of long jobs
  src/isotope.rs        formula parser and natural isotope patterns (M+0, M+1, …) from a built-in table
  src/ri.rs             retention index: RT ↔ RI calibration, marker search, RT shift, .calib file,
                        binary RICalculator record
  src/cmp.rs            saved components of the old version (.cmp): read and write
  src/search.rs         library search: spectrum and RI scores, one-to-one assignment of components
  src/batch.rs          batch quantification: many runs → substance × run table, quant ions
  src/batch/edit.rs     edits of the shown table (Edits, Act: drop or move runs and substances), hit index per run
  src/batch/export.rs   the export file set as contents in memory (results, groups, group statistics, settings, ions); md-app writes them
  src/normalize.rs      normalisation of the batch table (internal standard, total signal, pool, z-score)
  src/view.rs           view model of the tables: display rules (missing as 0, NaN sort, digits, axis index), batch view cached per generation; hit labels on the TIC, label overlap, LABELLED_PEAKS
  src/mid.rs            MID solver (corrections by sum formula, reference spectra, replicate mean)
  src/mid/report.rs     MidReport: the result of a MID scan, run and substance moves (shares `batch::pick` with the batch table), statistics, text exports (no egui)
  src/mid/              MID wizard pipeline (fragments, MID library, per-substance MIDs, group
                        statistics, the old text exports) and the FRAG "Export MIDs" of the batch table
  src/centroid.rs       centroided and nominal-mass spectra (NetCDF and MSL export)
  src/fastflight.rs     FastFlight TOF raw data (MsData.FFT/.FFC) import; src/tof.rs: its m/z calibration (.mzcalib)
  src/quadfit.rs        clean-room parabola fit that reproduces the old rounding of precise apexes (gsl-clean-room.md)
  src/legacy_numeric.rs the deliberate 3.3 roundings: C `(int)`, SG truncation, `roundSignificantFigures`, glibc pow(x, 2) and hypot
  src/library_add.rs    add a found component to a library, copy RTs into a target library
  src/spectrum.rs       spectrum normalisation (none, base peak, unit, relative)
  src/stats.rs          two-group biomarker table (Welch/paired t-test), five-number summary, PCA
  src/csv.rs            CSV writer for plotted series
  examples/perf.rs      performance smoke test on synthetic large runs
  examples/batch.rs     batch quantification timing on real files
  tests/                integration tests (example data, golden test against the old version)
crates/md-app/          egui desktop app (binary `metabolite-detector`), drawing only
  src/main.rs           entry point, `App` state and one call per UI area in `ui`
  src/menu.rs           top menu bar (File, Analysis, View, Help)
  src/status.rs         bottom status bar
  src/run_list.rs       left area: the open runs, the component list below them
  src/results.rs        the component list: sortable, works without a library, ↑/↓, CSV and .cmp (rows, order and filter in md-core)
  src/central.rs        centre: TIC/SIC plots, scan cursor, scan/RT entry, spectrum
  src/detected_peaks.rs Detected Peaks: ion traces of the selected component beside the spectrum
  src/view.rs           chromatogram view options (derivative, time unit, auto Y, score filter) and the run info window
  src/shortcuts.rs      keyboard shortcuts; menu items and keys run the same `Action`
  src/export_all.rs     File → Export all runs as NetCDF
  src/plot_export.rs    save a plot as PNG, SVG or PDF, or copy it; src/vector.rs writes the SVG/PDF from the frame's shapes
  src/tof_calibration.rs TOF calibration window (FastFlight m/z, saves .mzcalib)
  src/loading.rs        opening files, load polling, SIC jobs
  src/jobs.rs           `Job<T>`: a worker thread with `Progress`, `Cancel` and panic report; dropping it cancels. `Slot<T>` holds one job of a window (start, poll with repaint, cancel); new windows take a `Slot`
  src/calibration.rs    RI calibration window (alkane search, .calib load/save)
  src/isotopes.rs       natural isotope calculator window
  src/mid.rs            MID window: labelled/unlabelled runs, correction, result tree, bar plots, exports
  src/export.rs         spectrum normalisation picker, CSV and PNG export of the shown plots
  src/layout.rs         list/details split, main-window areas (`Panels`, `dock`), where tool windows first open; a test opens every
                        window and checks no control covers another, is cut off or reaches past
                        the screen
  src/identify.rs       library loading, hit labels on the TIC and the hit list of a component
  src/library_editor.rs library editor window (compounds, RI, quant ions, FRAG, comments; saves .mdlib)
  src/library_editor/add.rs  its "Add component" and "Copy RTs" dialogs (#80)
  src/batch.rs          batch quantification window and result table
  src/statistics.rs     biomarker statistics window: two-group table, box plots, PCA plots
  src/style.rs          app-wide look egui lacks (symbol font fallback, solid scroll bars, window margin); rules in DESIGN.md
  src/ui_tests.rs       glyph check for UI strings; ignored headless screenshot test (docs/ui-quality.md)
crates/md-cli/          headless batch quantification (binary `md-cli`), depends on `md-core` only
  src/main.rs           options, settings/library/calibration loading, CSV output
  tests/                runs the binary on the 15 example runs
docs/reverse-engineering/  evidence from the old binary (features, algorithms, file formats)
```

All analysis code lives in `md-core` so the GUI, CLI and tests share it.

Each crate follows the Cargo layout: `src/` for code, `tests/` for integration tests, `examples/` and
`benches/` for runnable measurements. A new crate goes under `crates/` and takes its versions, edition,
dependencies and lints from the workspace (`*.workspace = true`, `[lints] workspace = true`).

The crate boundary is the enforced architecture: `md-core` has no GUI dependency, so analysis code cannot
reach into the UI. SonarQube Cloud's Architecture view does not support Rust (only Java, JavaScript,
TypeScript, Python and C#), so it shows no architecture for this project.

### Code analysis

The Ubuntu job of `.github/workflows/ci.yml` runs SonarQube Cloud on `main` and on pull requests, reusing its own build: it measures coverage with
`cargo llvm-cov` (LCOV) and hands Sonar its Clippy JSON report; Sonar does not run Clippy itself. `sonar-project.properties` marks each crate's
`src/` as code and `tests/`, `examples/` and `benches/` as test code. Coverage counts `md-core` only: `md-app` just draws, so logic that
needs tests belongs in `md-core`. Automatic Analysis cannot analyze Rust,
so it must be off; the scan step is skipped until the repository secret `SONAR_TOKEN` exists.

## Performance rules

Large studies (hundreds of runs, tens of millions of points each) are the normal case.

1. **Compact storage.** Spectra are flat `f32` arrays with scan offsets (≈ 8 bytes per data point),
   matching the old `GCMSDiskScan<int, float>`. Accumulate sums in `f64`.
2. **Never block the UI.** File loading and every analysis step run on worker threads; the UI polls
   results and shows progress. Long jobs must be cancellable: every background job is a `jobs::Job`
   with the `Cancel` and `Progress` of `md_core::job`, which the command line uses too. Each operation has one function that takes the `Cancel` (callers without a job pass `&Cancel::default()`). Cancel is cooperative: a job checks between steps, and the long md-core steps take the `Cancel` themselves: `fastflight::bin_spectrum` (per scan), `andi::write_with_tic` (per variable), `Chromatogram::sic` / `ri::marker_peaks` (per scan) and `map_compounds_to_library` (per component); a plain file read still runs to its end. Opening runs can be cancelled from the status bar. Each job has a test that a set `Cancel` gives no result (for `write_with_tic`, `bin_spectrum`, `marker_peaks`, `sic` and `map_compounds_to_library` it is set inside the loop with `Cancel::after(n)`, md-core tests only); closing a job's window or dropping the `Job` discards its result. A cancel-checked parallel loop writes into a preallocated `Vec` instead of collecting `Option`s (the latter cost about 1.4 ms per SIC call).
3. **Draw proportional to pixels, not data.** Line plots go through `md_core::decimate`
   (min/max per pixel bucket of the visible range), so peaks stay visible at any zoom.
   The first derivative of a trace is cached per run (`Run::d_tic`, `Run::d_sic`), not computed per frame; both are freed when the derivative view is switched off.
4. **Parallel by default.** Per-scan, per-ion and per-file work uses rayon.
5. **Index once, query many.** m/z is sorted within each scan (binary search for SICs). Peak detection
   and deconvolution share one ion-major index per run (`IonTraces`) instead of rescanning spectra per ion.
   Per-point callbacks in these loops are generic (`impl FnMut`), not `&mut dyn FnMut`: the indirect call per point cost about 2 ms per index build (#396).
6. **Measure.** `examples/perf.rs` reports timings; each new algorithm adds a case, and PRs that
   touch hot paths report before/after numbers against the table below.
7. Memory-mapped or streamed access (like the old disk scan) is added when real runs exceed RAM;
   the cut-off is recorded as a `shortcut:` comment in `andi.rs`.

Release build, 36 000 scans × 500 ions (18 M points, 137 MiB). Median of 5 alternating runs of
`cargo run --release -p md-core --example perf` (the same `perf.rs` for all three builds), one session on
2026-10-22 on a noisy machine (single runs vary by about ±10 %). `f16f0a7` is 0.2.0, `93cbcb9` is main before
#395, `#395` is this fix. The "SIC median" row is the median of 15 calls; the single call also pays for touching the data.

| Operation | `f16f0a7` (0.2.0) | `93cbcb9` (before #395) | #395 |
|-----------|-------------------|-------------------------|------|
| TIC | 1.9 ms | 2.0 ms | 1.9 ms |
| SIC (one m/z), single call | 0.8 ms | 3.1 ms | 0.8 ms |
| SIC (one m/z), median of 15 calls | 0.20 ms | 1.56 ms | 0.17 ms |
| Decimate TIC for drawing | 0.04 ms | 0.05 ms | 0.04 ms |
| Overlay frame: decimate TIC + SIC of 15 runs | 0.66 ms | 0.83 ms | 0.63 ms |
| Ion traces (index, 500 traces) | 19 ms | 21 ms | 22 ms |
| Peaks in all 500 traces (300 000 peaks) | 17 ms | 15 ms | 15 ms |
| Deconvolution (12 000 substances of 25 ions) | 19 ms | 20 ms | 20 ms |
| Library search, 1 000 components × 10 000 entries, cutoff 0.8 | 33 ms | 33 ms | 32 ms |
| Library search, same, cutoff 0 (9 M hits kept, like the old GUI) | 243 ms | 70 ms (#370) | 69 ms |
| `run_study`, 1 run targeted, rescan on (file read, traces, peaks, deconvolution, search) | - | - | 419 ms (#405, see below) |
| `run_study`, 1 run untargeted | - | - | 441 ms (#405) |
| `run_study`, 4 runs targeted, 4 in parallel | - | - | 1 277 ms (#405) |
| Exact-mass run, 20 000 m/z (sparse index ≈ 0.15 GiB): traces / peaks + reorder / deconvolution | not measured (no case in `f16f0a7`; dense peak memory 2.8 GiB) | 459 / 461 / 48 ms (#369, earlier session) | not re-measured |

#405 (one session, 6 alternating pairs of `cargo run --release -p md-core --example perf`, medians, `fab16d3` (main) before and this branch after; the batch rows came from a `scan_chromatogram` that is test-only now, so they are replaced by the `run_study` rows). The machine was noisy: the TIC row moved although its code did not change, so differences below about 10 % are not significant.

| Operation | main (before) | #405 (after) |
|---|---|---|
| `run_study`, 1 run targeted, rescan on | 426 ms | 419 ms |
| `run_study`, 1 run untargeted | 437 ms | 441 ms |
| `run_study`, 4 runs targeted | 1 267 ms | 1 277 ms |
| Ion traces (500) | 26.0 ms | 26.7 ms |
| Peaks in all traces | 17.2 ms | 18.7 ms |
| Deconvolution | 23.1 ms | 23.8 ms |

The removed `fresh.clone()` copied one peak list per run; that is not measurable at this size. The gain is the cleaner interface and no copy on nominal traces.

SIC was the regression: `sic_cancellable` (now `sic`, #378, `5f7cf5f`) collected `Option`s in its parallel loop, which cost about
1.4 ms per call. The fix keeps the per-scan `Cancel` check, writes into a preallocated `Vec` and returns `None` at the end;
`map_compounds_to_library` now collects its hit lists directly instead of collecting `Option` values, which takes about 8 % off the batch run.
Ion traces (#396): the index build called its per-point callback through `dyn FnMut`, which the dense store could not
inline (#369). The callback is generic now. Median of 7 alternating runs in one session, 36 000 × 500: ion traces 23.0 ms
before and 20.9 ms after (`f16f0a7`: 20.7 ms); the full batch run (`scan_chromatogram`, median of 5) 102.8 ms before and
99.3 ms after (`f16f0a7`: 95.1 ms). `examples/perf.rs`, median of 7: ion traces 22.0 ms before and 20.3 ms after
(`f16f0a7`: 23.9 ms); one batch run 315 ms before and 314 ms after (`f16f0a7`: 303 ms; single runs vary by about ±10 %).
The step `4f136f9` (#377) does not show in this session: 98.6 ms at `d993ca5` against 97.4 ms there (median of 7).
Of the batch gap (4.2 ms, about 4 %), the assignment heap of #370 in `map_compounds_to_library` is a partial cause, about 1.9 ms:
the search over the 12 000 components of the run takes 21.2 ms here against 19.3 ms in `f16f0a7` (3 runs). The heap stays on purpose: the pair
list it replaced needed 16 bytes per hit (about 144 MB at cutoff 0). The rest of the gap, about 2 ms, is not explained.
Tracked in #396.

Real runs (TargetSearchData, 4 400 scans, ~745 000 points, 6 MB `.cdf`): read + TIC 13–25 ms per
file; all 15 runs in 0.27 s; ion traces + peak detection 4–8 ms per file, 630–950 peaks in 236 traces;
with deconvolution 5–10 ms per file, 13–20 components of 15–21 ions
(`cargo test --release -p md-core --test example_data -- --nocapture`).
The app opens all 15 as an overlay in 0.06–0.21 s (one loader thread per file; status bar shows the time).
The whole file is read into memory first because `netcdf3` seeks per record (1.8 s per file on a
plain `File`).

## Peak detection

`md_core::peaks` ports the old `DefaultPeakDetector` ([algorithms.md](reverse-engineering/algorithms.md)
§2–3): `IonTraces::new` builds the ion-major index once per run, `IonTraces::detect_peaks` runs the
detector on every trace in parallel. The index is dense (masses × scans × 4 bytes; the example runs:
about 0.1 GiB) or, when the run fills at most 1 in 4 cells (`SPARSE_FILL_RATIO`, typical for exact-mass
runs), sparse: per row the scans with a point (8 bytes each), every other scan reads the missing value.
Readers call `IonTraces::trace(row, scans, buf)`, which returns the dense slice or fills `buf` (one per
rayon worker; peak detection keeps it thread-local because `map_init` would allocate it per row), and
`IonTraces::value(row, scan)` for single points. The ratio is measured with `examples/perf.rs` on
36 000 scans × 500 points per scan spread over 1, 2, 4, 8 and 40 exact m/z per mass: at 1 in 2 sparse
saves nothing, at 1 in 4 it saves a third of the peak memory for about a quarter more time, at 1 in 40
(20 000 traces) the peak memory falls from 2.8 GiB to 0.3 GiB at the same detection time.
`reorder_peaks` keeps a running f64 sum per profile, so a join adds one member instead of re-adding all.
`IonTraces::new` returns `TooLarge` above `MAX_TRACE_BYTES` (8 GiB) before allocating, so a
corrupt file cannot exhaust memory, and above `MAX_TRACE_CELLS` (8 Gi masses × scans, sparse too: detection
reads every cell), so a corrupt exact-mass file with a distinct m/z per point cannot run for minutes. The app shows that error in place of the run. `PeakSettings::default()` holds the old `AN_*` defaults.
`matched_filter` ports the old filter, which the deconvolution also uses (§4.5).

`PeakSettings::detector` selects the old "Matched Filtering" method instead (`Detector::Matched`, old
`MatchedPeakDetector`, §3.3). The old app offers it only in its *Calibration Setter* window, with peak
width 10 and threshold 3 as the values after switching (`Detector::MATCHED`); the main window always uses
the default detector, and the old settings file has no field for it. Here it is a field of `GcmsSettings`
(not saved in the old file) and a choice in the app's settings window, so the app, the batch scan and the
deconvolution all use the chosen detector. In the old code the matched peaks hold the filtered values and
the deconvolution reads those, so `IonTraces::detect_peaks` overwrites each trace in place with its filtered
values (one row copy per worker, no second index). The filter is computed per trace in parallel with one
precomputed kernel; `detect_peaks_matched` runs the same detector on a single trace (RI markers).

`md_core::analysis` puts the chain together once: `analyze_run` (index, peaks, deconvolution) for the app
and the batch, and its peak stage `find_run_peaks` for the rescan, which keeps the index and deconvolutes
only small windows. `PeakPass` names the one difference: `Reordered` is the main window and batch
(`detect_all_peaks`, old `detectAllPeaks(reorder=true)`), `Fresh` is the rescan (`detect_peaks` only,
like the old `ResultSeries::scanIons`, algorithms.md §4.0). Both passes also return the `Fresh` peak list
(`RunPeaks::fresh`, `RunAnalysis::fresh`), which the batch keeps for its rescan. `cancel` is looked at
between the stages (index, peaks, deconvolution; the batch library search then checks before each component), so a
cancelled run stops at the end of its current stage. The batch passes its run by value so the
spectra are dropped once the index is built; the app lends its run and keeps it for the plots.

![Detected peaks on four traces of 7235eg04.cdf](images/peaks-7235eg04.svg)

Conscious deviations from the old app; the golden test (below) shows that none of them changes the
peaks of the example runs:

- **Traces per nominal mass.** Not a deviation: `IonTraces::new` centroids every scan like the old CDF
  import (§4.1), then sums the centroids per nominal mass. On unit-resolution runs with integer
  intensities centroiding changes nothing. Runs with `Chromatogram::exact_masses` (FastFlight, caches
  with DataType 0 or 1) get exact-mass traces and `reorder_peaks` (§4.2) through
  `IonTraces::detect_all_peaks`, like the old app. Components keep the exact masses in
  `exact_spectrum`; `spectrum`, `.cmp` and the library search still use nominal masses (shortcut in
  `deconv.rs`).
- **Noise factor 1.0.** `GCMSAnalyzer` sets the detector's noise factor from
  `DataProvider::getGlobalNoiseFactor` (§2.3). On all 15 example runs the old app logs
  "Error determining nf. Using nf=1", so 1.0 is what it uses there. Data where the old app finds another
  noise factor would differ, e.g. FastFlight caches with DataType 1, whose per-scan
  `analyzeBaseline(values, 100)` statistics set the noise factor (§2.3). The noise factor and MAD of
  `analyzeBaseline` are not computed because the detector never reads them.
- **f32 intensities instead of int.** Only the SG smoothing and the matched filter output are truncated
  to integers like the old int filters (`legacy_numeric::sg_trunc` and `c_int`, the x86 `(int)`
  with `INT_MIN` out of range). S/N and sharpness are not rounded to 8 significant digits.
- **Two `analyzeBaseline`.** The old program has two functions of that name: the peak detector's
  (`peaks::trace_baseline`, segments of 15 scans accepted by slope) and the spectrum one
  (`centroid::spectrum_baseline`, segments of 25 points accepted by mean crossings). They stay
  separate; a test shows a ramp where they disagree.
- **Precise apex: not a deviation.** The old app fits the parabola with `gsl_multifit_linear` in absolute
  bin coordinates; at an exact half-bin tie and on a flat top (three equal values) the fit's rounding
  decides the result. `quadfit` rebuilds that fit bit for bit from black-box outputs of the old binary,
  without GSL code (`docs/reverse-engineering/gsl-clean-room.md`), so ties and flat tops land where the
  old app puts them, for every bin count while bins × scan stays below 9.4e7. Above, the x² entry may
  round differently from the old `pow` (shortcut in `legacy_numeric::glibc_pow2`).
- **No progress or cancel per 100 traces.** A whole run takes milliseconds (table above).
- **Border extension off** like the old default (`setPeakBorderExtensionEnabled` has no known caller),
  so the extension code is not ported.

## Deconvolution

`md_core::deconv::deconvolute` ports the old `deconvoluteChromatogram` ([algorithms.md](reverse-engineering/algorithms.md)
§4.4–4.9): the apex intensities of all peaks go into a histogram of `AN_NUMBER_BINS` bins per scan, the
matched filter (`AN_DECONVOLUTION_WIDTH` scans) marks regions above 5, and each region becomes a component.
Its model peak is the sum of the sharpest clean peaks (quality ≤ 10, at most a quarter); every ion trace is
fitted as `c0 + c1 * scan + c2 * model`, and `c2 * max(model)` is that ion's height in the spectrum.
Components are fitted in parallel. `DeconvSettings::default()` holds the old `AN_*` defaults; the golden tests
use it. Deviation (#279): the defaults of app and CLI (`GcmsSettings::default()`) set the width to automatic,
`deconv::auto_width`: 0.5 s at the run's median scan interval, at least the old 1 scan, at most 100 and a tenth of the run (#312). 3.3 always uses
1 scan, which cuts the substances of a 20 scans/s run into pieces (#242). A width from a loaded settings file, the
imported 3.3 settings or the saved app settings (`AN_DECONVOLUTION_WIDTH=auto` or a number) is a set one and wins;
the old file format has no "auto", so **Save…** writes the width of the shown run (1 without one). The app runs
peak detection and deconvolution after loading and marks each component on the TIC; a click near a dot
shows its deconvoluted spectrum.

The old quirks are kept on purpose; with them, all 246 old components match exactly (golden test below): the first
and last bin of a region are added twice, so a peak can appear twice in a component and counts twice
towards `AN_MIN_PEAK_NUMBER`; the model peak sums raw intensities including the baseline, so all heights
of a spectrum share one scale factor.

Conscious deviations:

- **No precise apex → bin 0.** Not a deviation: where the parabola fit was refused (apex at the peak
  edge, vertex a scan or more away), the old `interpolatedMaxScan` keeps its initial 0, so the histogram
  and `reorderPeaks` put the peak into bin 0 (`Peak::interpolated_max_scan`). Seen with the matched
  detector, whose peaks can start at their apex.
- **Collinear fits drop the peak.** The old code uses an SVD (`gsl_multifit_linear`), which also returns a
  minimum-norm fit when model and scan are collinear; here those peaks are skipped.
- **Filter edges are zero-padded** (old edge handling unknown, §4.11).
- **Not a deviation: `roundSignificantFigures`.** `legacy_numeric::round_sig8` follows the old code
  step by step: round half up to one digit more, then floor, so a share below 1 keeps 8 significant
  digits and one from 1 up keeps 9. With it every model share, height and integral of the stage
  fixtures matches (a plain round to 8 digits missed 1938 of 4366 shares by one unit).
- **No `reorderPeaks`** (unit-mass data only, see peak detection).
- `AN_MIN_HEIGHT_BP` has no static value in the old binary (zero-initialised), so the default keeps all ions.

## Retention index

`md_core::ri::RiCalibration` ports the old `RICalculator` ([algorithms.md](reverse-engineering/algorithms.md) §6):
RT → RI is piecewise linear through the calibration points, and outside them the first or last segment
continues (§6.2). `rt` is the inverse (§6.3), `component_ris` gives the RI of each component at its apex
scan (§6.4). Times are seconds like `Chromatogram`; the `.calib` file holds milliseconds.

Calibration from a run: `marker_peaks` runs the matched detector (§3.3) on the marker ion trace (± 0.5)
and takes the highest raw scan of each peak. `series_points` gives the n highest peaks, in time order, the
RIs of an alkane series (C10, C11, … → 1000, 1100, …); `window_points` takes the highest peak in a known
time window per marker. The app's *Analysis → RI calibration* window finds alkanes in the active run
(m/z 57, C10–C36 by default) or loads/saves a `.calib`; the calibration applies to all open runs and to the
library search, and the status bar shows the RI of the scan or selected component.

Calibration mix (old wizard, [algorithms.md](reverse-engineering/algorithms.md) §6.8): the user picks library
compounds with RI; the calibration run's components are searched against them by spectrum only, and
`match_mix` assigns them as the old `matchCompounds` does (mix by ascending RI, each takes the next component
in time order that has it as a hit). The result is an editable table (use / RT / RI); substances not found
are unused rows. A change takes effect when the field is left or the drag ends and builds a new
`RiCalibration` (invalid tables keep the old one and show the error). Each run can take over the table as
its own calibration and save it as `.calib`; runs without one use the table. Deviations: the old table offered
the components' times in a combo box, here the time is typed; substances with the same RI are both kept.
*RT only* (old "Use Retention time for chromatogram alignment (not recommended)", with its warning):
`RiCalibration::identity()` (RI = time in s) and `PreparedLibrary::new_rt_only` (library RI = library RT in s)
make the library search and the reference shift compare times. The RT/1000 fallback of §5.6 is still not
ported for the normal mode.

Checked against TargetSearch on the example data (`tests/example_data.rs`): its FAME markers (m/z 87 in
the windows of `rimLimits.txt`) give the RI of every scan in `RI_*.txt` to within 0.01 in 14 of 15 runs.
In `7235eg06` TargetSearch took a peak at 321.1 s as the middle marker, where m/z 87 reaches 36 512;
the marker peak at 311.1 s, as in all other runs, reaches 221 839. Filter width 40 scans, threshold 100.

RT shift correction (§6.5): `pick_reference` chooses the substance with the best library score for the
reference near its RI, `shifted_to` shifts the calibration so it lands on the library RI.
`md_core::search::shift_to_reference` combines them with the library search: window `3 · LS_RI_DIFF`, score
= overall score of the reference among the component's hits. Because that score contains the RI score, a run
shifted by more than about half a sigma (10 RI) no longer reaches the default cutoff 0.8; the old code
behaves the same. In the app the calibration window offers a reference substance (library compounds with RI)
and the cutoff once a calibration and a library are loaded. Each run's identification job then shifts the
calibration onto the reference, searches with the shifted RIs and keeps the shifted calibration for the
status bar; the window lists the shift per run, or that the reference was not found (then the run uses the
unshifted calibration). Reference and cutoff are part of the `run.ident` key. Tests: a synthetic run 10 s
late (`search.rs`) and the example runs made 4 s late, where Glycine 3TMS lands on its `lib_aa` RI again and
the shift grows by exactly 4 s (`tests/example_data.rs`).

One precedence rule for the RI a run shows (`md_core::ri::applicable`, called by `Calibration::for_run`; the
component list, the calibration window and identify labels read its result): the identification's calibration
(shifted, else the one it was searched with) applies, also while a new identification runs; without one, the
base applies. RT-only: without an identification nothing; with one, its calibration maps RT to the shifted
time. Table-driven test in `ri.rs`. Behaviour change (#407): in RT-only mode the detail panel and the run list of identified runs now show the
time as RI, where before they showed nothing; the component list is unchanged.

Conscious deviations:

- **At least two points, RI rising with RT.** The old code accepted any table and returned -1 below two
  points; a falling table has no inverse. A slope outside 1e-9 to 1e9 RI per second would overflow on
  extrapolation. Such calibrations and `.calib` files are rejected with an error.
- **Above the last point** the old code dereferenced `end()` (undefined, §6.2); here the last segment
  continues, as for every other time above the range.
- **The shift does not stack.** The old code subtracted the previous shift once more (§6.5), so only
  the first correction was exact.
- **`.calib`** is read as `<rt in ms> <RI>` per line (any whitespace, blank lines skipped; the ms unit
  is inferred, file-formats.md §4.1) and written with full precision instead of 6 significant digits.
  Reading is a trust boundary: at most 1 MiB, UTF-8, two finite numbers per line, errors name the line.
- **Marker apex** is the highest raw scan of the matched-filter peak; the old dialog's choice is unknown (§6.7).

## Library search

`md_core::search` ports the old `LibrarySearch` ([algorithms.md](reverse-engineering/algorithms.md) §5).
`PreparedLibrary::new` bins, mass-filters and scales every library spectrum once (mass² · intensity^0.6 with the
default `LS_SCALE_FACTOR_*`), like the old pre-scaled cache. `library_hits` scores one spectrum against every
entry: the pure score is a cosine over all peaks of both spectra, the impure score over the library peaks only
(peaks without a partner are weighted with `LS_SCALE_FACTOR_NO_MATCH`), the spectrum score is
`0.5 · pure + 0.5 · impure`, and with an RI the overall score is `(spec² · ri)^(1/3)`, where the RI score is a
Gaussian with sigma `LS_RI_DIFF` that is zero beyond 3 sigma. `map_compounds_to_library` searches all
components in parallel (`map_compounds_to_library` checks `Cancel` before each component and returns `None` once set; the batch scan uses it) and then assigns library compounds greedily by best overall score, each library compound
and each component at most once. The assignment does not build or sort the list of all (component, hit) pairs:
each component's hits are already best first, so a heap over the components' current best pairs takes them in
the same order (ties: lower library index, then lower component index) and gives the same result bit for bit
(`assign_best_first`, checked against the sorted pair list in a test). All hits stay, the colouring needs them.
`SearchSettings::default()` holds the old `LS_*` defaults.
With nominal masses each pair is one pass over the library peaks with a dense lookup into the component
spectrum, so no fragment pre-filter is needed for speed (table above). For NIST text libraries (`.sdf`,
`library::is_nist`) the app and `md-cli` turn on the old NIST pre-filter
(`PreparedLibrary::with_fragment_filter`, file-formats.md §6.3): only entries that share more than
`min_equal_fragments` (3) of the `max_fragments` (5) strongest fragments of the component are scored, as in
the old NIST search. Both are in the settings dialog ("Number of fragments", "Minimum number of identical
frag."; old `SQL_*`).

File → Open NIST library loads a second library next to the normal one, always with the pre-filter (old
`SQLLibrary`). The component panel then has "NIST search": the selected component is searched in it on a
worker thread (about 10 ms for 300 000 entries with the default pre-filter; a low minimum scores almost all
entries; the button waits meanwhile) and its NIST hits are listed below the library hits, which stay
as they are (old "NIST search..."). The NIST hits belong to that component and spectrum; a re-analysis hides
them.

The app opens a library with File → Open library (or a `.mdlib`/`.lbr`/`.msl`/`.msp`/`.sdf` path on the command line), identifies every
open run on worker threads, writes the assigned name at each component on the TIC (green above 0.9, orange
above 0.7, red otherwise, like the old app) and lists the best hits of the selected component with overall,
spectrum and RI score. With an RI calibration the component RIs (`component_ris`) go into the search; setting,
loading or clearing the calibration re-identifies every run, because `run.ident` remembers the library and
calibration (`Arc`s, compared by pointer) it was computed with. `tests/example_data.rs` checks this with
`lib_aa.lbr` (`MD_LBR_SAMPLE`): the example runs have no alkanes, so their alkane calibration is bridged over
Valine 2TMS and Threonine 3TMS (FAME RI from TargetSearch's `library.txt` → time in the run); Leucine 2TMS,
Glycine 3TMS and Serine 3TMS are then assigned with RI scores of 0.94–0.99.

Conscious deviations:

- **No RI → spectrum only.** Without a calibration, or when the library entry has no RI, the overall score is
  the spectrum score alone (with `use_spec` off there is then no hit). The old code then compares retention
  times / 1000 with the same sigma instead (§5.6), and with `LS_USE_RI` on it drops every entry more than 60 s
  away. That is not ported: times of different instruments and temperature programs are not comparable. The
  sample library `lib_aa.lbr` lists Glycine 3TMS at 418 s, the example runs elute it at 313 s, so the old rule
  would remove the correct hit there. The old NIST search drops entries without RI and finds nothing for a
  component without RI; here those are scored by spectrum as well.
- **Nominal mass.** Library spectra are rounded to nominal mass and intensities of one mass are summed, because
  the old `getBinnedCentroidCompound` is not decoded (§5.9) and components are unit-mass. Matching is then exact,
  which makes the old 0.01 m/z tolerance and its double-match quirk (§5.3) irrelevant. Negative intensities
  and masses below 0.5 are dropped.
- **Hits need a score above zero.** The old search keeps pairs with score 0 and may even assign them; here they are
  no hits. The cutoff is applied during the search; the app searches with 0 like the old GUI. A higher cutoff
  yields the same assignments above it, because the greedy assignment takes pairs in descending score.
- **Ties are deterministic**: lower library index, then lower component index (old: pointer order).
- The Wave score variants are not ported (§5.9). `LS_USE_SCALED_LIB` is only stored: 3.3 reads it nowhere either (§4b.4).
  `LS_IDENTIFICATION_ALGORITHM` is not ported either: the old app never writes it, so only algorithm 0 (the
  ported one) runs there (file-formats.md §5).
- **Settings are not kept between sessions.** The old app stores the NIST path (`rc_lib_nist_file`) and the
  `SQL_*` values with QSettings; the rewrite keeps no app preferences yet, so the NIST library is opened per
  session and the two values hold until the app closes (loading an old settings file keeps them).

## Batch quantification

`md_core::batch` finds the library substances in many runs (old `BatchQuantificator`, F10–F12) in two
passes behind one call, `run_study(&Study, progress, cancel) -> StudyResult`, which the batch window, the MID
window and `md-cli` all use (`Study`: runs, library, scan settings (which carry the quant settings), rescan switch, runs in memory,
default `DEFAULT_PARALLEL` = 4). Pass 1 reads each file, runs ion traces, peaks, deconvolution and library search (with RI if
a calibration is given) and keeps only the assigned substances with their ion integrals. Files go through
`parallel` worker threads, so at most that many runs are in memory; each run is analysed with rayon as
well. Progress and cancel are `md_core::job::Progress` and `Cancel`; cancelling stops the runs in work at their
next stage and starts no further file; unreadable files are listed and skipped. Each `RunResult` carries its
`input` (position in `Study::runs`), so callers match results to inputs by position, not by path. The quant ions
of the study are kept in `StudyResult::quant`; `table(&StudyResult, library, TableSettings, quant)` builds the substance
× run table and is cheap, so the app re-filters without a rescan. `quant` (`Option<&QuantSettings>`, `None` = the scan's
own) picks other quant ions after the scan, but only ions the scan already has: a scan with a rescan keeps its table
until it is scanned again (the app says so; #411):

- **Quant ions** per substance: the library's manual ions, else its automatic ones, else the old
  `autoQuantificationIonAssignment` (§7.2, fallback quirk kept) on the best-scoring run. The value of a
  cell is the summed integral of the quant ions found in that run.
- **Replicate groups**: the first capture group (else the whole match) of a regex on the file name.
  A sample table (`SampleTable`, CSV/TSV; file name → group; extension and case do not matter) wins for the
  runs it names; the CLI takes it with `--sample-table`, `--sample-file-column` and `--sample-group-column`.
- **Filters**: minimum quant ions found per cell, minimum reproducibility (share of runs in a group),
  minimum S/N of an ion and maximum discrepancy index (mean peak quality of the ions of a cell). A
  substance without any value left drops out of the table.
- **Untargeted** (`ScanSettings::untargeted`, then `add_unknowns`): components without library hit are kept
  per run and pooled over all runs into a temporary library (`Unknown n`, the mean spectrum of its members);
  `table` takes it as `extra`. The components of a run are scored in parallel against the library as it
  was before the run; 12 000 unknowns in 20 runs take 0.07 s on 32 logical cores and 0.47 s on one thread (`examples/perf.rs`).
- **Rescan** (`rescan_runs`): a substance found in at least a quarter of the runs gets its missing quant ions
  from the ion traces of the runs' files, within ΔRI (`SearchSettings::ri_diff`) around the expected
  position (old `ResultSeries::scanIons`). It reuses pass 1's `Fresh` peak list of each run (#363; about
  750 peaks or 27 KiB per example run), not its ion index: a run whose list has no peak of a missing ion
  in its window is not read again, and a run read again skips the peak detection (the matched detector
  rewrites the traces it detects on, so it detects again). The result is the same as with a full second
  detection (`study_equals_the_passes_with_a_full_rescan` on the 15 example runs). The local
  deconvolution of each missing ion is most of the rescan time; caching the ion index is not done.
- **RT shift** (`ScanSettings::shift`): the calibration is shifted per run onto a reference substance.

What the old binary does and what is own choice: `docs/reverse-engineering/batch-quantification.md`.

**Result view** (#83): the window can remove runs and substances, move runs, show the retention time of
each hit (`Row::rts_s`), the total intensity per run (`Table::total_intensities`) and per replicate group
the old statistics columns (`stats::summary`, `stats::student`). `Table::select_runs` and
`Table::without_rows` do the edits; runs are removed before the normalisation, substances after it, so
the CSV export and the statistics window see the edited table. Old formulas (from the 3.3 binary):
"Mean" and "Rel. Std. Err." (sample SD / sqrt(n) / mean, a fraction) use only values that are not 0, NaN or
infinite; "#" counts values above 0; "P (blank)" is Student's t-test with pooled variance, two-sided, with
missing values as 0, against the group named `BLANK`.
The old table gets the internal standard and summed-signal normalisation (`normalizationPerformed` →
`BatchQuantificationStatisticsTableModel::setNormalizationSeries` → `setNormalizationBySum` /
`setNormalizationSeries`), never the z-score, which is only a flag of `getTableData`. Its cells are
`QVariant(double)`, which Qt 5.15 shows with 6 significant digits. So the rewrite shows the mean with 6
significant digits (`significant`, without exponent) and offers no "Group statistics" while the z-score is
on (#189). Not checked: whether the old pool normalisation (`ResultSeries::setNormalizer`) reaches the group
statistics; the rewrite applies it.

**Details behind a cell** (#160): `Found::alternatives` keeps the other library hits of the component
(`md_core::batch::hit_list`, assigned hit first); the Values view shows them as a tooltip, the assigned one in
red, like the old tree. "Quant ion lines" adds a line per quant ion under each substance with its raw
integral (`ion_integrals`, the same rule as the cell: largest peak with at least `MIN_SN`). An export also
writes `FILE.settings.csv` (`Setting,Value`: the batch filters and `QUANT_*`, the normalisation and the number
of removed runs and substances; `TableSettings::write_csv`) and, on request ("Export quant ions" /
`md-cli --ion-integrals`), `FILE.ions.csv` (`Substance,Quant_ion,<run>...`, raw integrals, empty = not found;
`Table::write_ions_csv`); the settings file also names the library, RT-only mode, search cut-off and the
number of scanned runs, unknown substances and ions added by the rescan (md-cli: library, settings and
calibration file, shift reference, `--untargeted`, `--rescan`). `scripts/r/load_results.R` loads both next to the
table as `settings` and `ions`. The scan result keeps an index per run (`Scan::index`), so the tooltip and the ion
lines look up a hit in constant time. The md-cli file names files, not paths. The last column of "Group statistics" is the one-way ANOVA p-value over all groups (`stats::anova`, #203; formula in `reverse-engineering/batch-quantification.md`).

All 15 example runs against the TargetSearchData library (12 entries, converted to MSL), measured
2026-10-20 on a loaded machine: 1.1 s with 1 run in memory (33 MiB peak working set), 0.17 s with 4
(90 MiB), 0.10 s with 15 (189 MiB). 150 runs (the 15 ten times) with 4 take 2.6 s and stay at 82 MiB, so
memory depends on the runs in memory, not on the study size (`examples/batch.rs`). Reading the file
(~27 ms per run) is the largest share. Untargeted, the 15 runs give 97 unknown substances besides the 5
library substances; with the rescan 380 of 1530 cells of 102 substances are filled (435 ions added, 2 ms
pooling, 1.1 s rescan with 4 in parallel), measured 2026-10-20 with `examples/batch.rs` after #110.
Reusing the pass-1 peaks (#363) takes the rescan of the 15 runs with 4 in parallel from 80 to 66 ms
(5 substances) and from 186 to 159 ms (untargeted, 474 ions added), medians of 15 in-process runs on a
loaded machine, 2026-10-22.
`examples/batch.rs` on all 15 runs, 4 in parallel, median of 5, same session as the table above (#360): 5 substances, no rescan,
117 ms before M8 (`f16f0a7`, scan) and 108 ms after (`2fb4769`, `run_study`); untargeted with rescan
375 ms before (scan + pooling + rescan, summed per invocation, then median) and 313 ms after (`run_study` does all three).

### Command line (`md-cli`)

`md-cli` runs the same two passes without a window, for big studies on servers or clusters:

```sh
cargo run --release -p md-cli -- --library lib.msl --settings algorithms.settings --calibration run.calib \
    --replicates '^(.+)_\d+\.cdf$' -j 8 -o table.csv study/*.cdf
```

Runs are `.cdf`/`.nc` or old `.bin` files (shell globs expand them). `--settings` is the old settings file
(`GcmsSettings::from_bytes`): its `AN_*`/`LS_*` values drive peak detection, deconvolution and search, its
`QUANT_*` values are the quant defaults, which `--ions`, `--min-dist` and `--quality-cutoff` override.
`--replicates`, `--min-reproducibility`, `--min-ions`, `--min-sn`, `--max-discrepancy` and `--exclude-ions`
(`QUANT_EXCL_IONS`) are the batch window's filters; `--untargeted`, `--rescan` and
`--shift-reference` (needs `-c`; `--shift-cutoff`) switch on the passes above. The CSV is
`Table::write_csv` (UTF-8 with byte order mark for Excel, `Substance,RT_s,Quant_ions,<run>...`, empty cell = not
quantified; names that start with `=`, `+`, `-` or `@` get a leading `'`); with replicate groups, `-o FILE`
and the app's export also write `FILE.groups.csv` (`Run,Group`, `Table::write_groups_file`) and
`FILE.groupstats.csv` (`Table::write_group_stats_file`: per group `Mean_`, `Rel_Std_Err_`, `Hits_` and with a group
`BLANK` (the app: the group named in "Blank group"; md-cli: `BLANK`) `P_blank_`, last `P_ANOVA`, unrounded, empty = not defined; the same numbers as the "Group statistics" view,
which takes them from `Table::group_stats`; the app writes it not for z-scores, #222; the old `exportData` writes the shown table cell by cell, so its group
view had the p-value column too, #222); `-o FILE` always
writes `FILE.settings.csv` and, with `--ion-integrals`, `FILE.ions.csv` (see "Details behind a cell"); progress and
failed runs go to stderr, and the exit code is 2 if the table is complete except for unreadable runs.
The default CSV format (`,` and decimal points) never changes, so scripts keep working; a spreadsheet with German
locale reads `.` as thousands mark (`293.808` becomes 293808, #249). `--excel-de` (the app: "Excel (German)")
rewrites every file of the export with `csv::to_excel_de`: `;` separators and decimal commas in the number
columns (also `RT_s`; `Quant_ions` is text with spaces and stays); header lines, name columns (so a group `0.5` stays)
and the values of `FILE.settings.csv` (text for `load_results.R`) keep their text, quoting as before. `load_results.R`
reads both; run with R 4.5.2 on the 15 example runs with groups (#322): `data` loads, groups `04` stay text, a leading
guard `'` is stripped, numbers are numeric (a test in `crates/md-cli/tests/example_data.rs` repeats this and skips
without `Rscript`). A group name like `01` is still quoted (`"01"`), but there is **no encoding that keeps it text on a plain
open**: Calc (de-DE, checked headless) reads `"01"`, `01` and `+01` as 1 and shows `'01` with the apostrophe; only `="01"` stays text, a
formula that the guard on `=` would have to let through (and `load_results.R` could not read back). Text stays only with "Quoted field as text" in
the import dialog, so the export dialog says so when a group like `01` exists with "Excel (German)" on (Excel behaves the same as far as
known, not checked); CRLF files keep CRLF. Checked in LibreOffice Calc headless (numbers stay numbers, names text).
Without `--settings` and `--calibration` the table equals the app's batch window with default settings
and no RI calibration. `crates/md-cli/tests/example_data.rs` checks this on the 15 example runs.

## Old result files

`md_core::cmp` reads and writes the old `.cmp` files (saved components, version 4 of MetaboliteDetector 3.3;
[file-formats.md](reverse-engineering/file-formats.md) §3) with every field, including the ion peaks and the
model peak, as `library::LibraryCompound` records (the reader is shared with `.lbr`). `LibraryCompound::component` turns one into the `deconv::Component` the app shows. A `.calib` is
read by `RiCalibration::read_calib`; the binary calibration the old RI wizard stores in a run's `.idx` by
`diskscan::read_calibration` (`ri::StoredCalibration`, values kept in ms so it writes back the same bytes).
Opening a run in the app uses `cmp::saved_components` and `ri::saved_calibration`: components from
`<run>.cmp` replace the analysis until the settings change, and `<run>.calib` (else the `.idx` calibration)
becomes the run's own calibration. A `.calib` has no RT shift, so it drops a shift stored in the `.idx`;
it wins because it is the file the user saved last (the app's "Save .calib" suggests that name). A `.cmp`
whose apex times do not match the run's scan times, a broken file and versions 1 to 3 of `.cmp` (no sample
to decode them from) only give a warning: the run opens with its own analysis. The ion peaks of the file
become the run's `Peak` list (`cmp::SavedComponents`), so the component panel shows S/N, sharpness and
discrepancy for loaded components too (#126); a saved peak has no ion-trace row (`Peak::row` is 0).

**Save and load (#107).** The component list (left of the main window, #251) writes the components of the active run
(`cmp::compounds_of_components`, `write_cmp`; default name `<run>.cmp` next to the run, where opening finds
it again) and loads a `.cmp` in place of them (`cmp::components_from_file`; selection and hits are reset and
the library search runs again). Like the old file, it keeps one ion peak per mass (the sharpest, as the panel
counts), ascending by mass, so a loaded component lists its peaks in that order. The model peak has the scan
range and values only and the name stays empty: the old program accepts that and gives the component a
numeric name when it reads it (Observed, `make-reference.sh --reread`; apart from that name, `time_stamp` and
the check value it rewrites our files byte for byte). Settings changes re-analyse and replace loaded components.
A `.cmp` is also a batch input (`batch::run_study`, file picker and `md-cli`): its components are searched
as saved, with RI and shift correction from the apex times; there is no raw data, so `rescan_runs` skips it.
A saved component keeps one peak per mass, so where a run has two peaks of a quant mass and only the second
passes `--min-sn`, the table of the `.cmp` lacks that value.
The batch also saves each run it deconvolutes as `.cmp` (`ScanSettings::save_cmp`; batch window "Save
components as .cmp", `md-cli --save-cmp DIR`), in place of the old Compound Detection Wizard
(`RICalibrationWizard` → `CompoundDetector::detect`, 0x3f4420: per run `detectAllPeaks`, deconvolution,
`calculateRI` with the run's calibration, then `<base>.cmp` next to the run; a failed run gives "Error
detecting compounds of file:" and the others go on). The file is written in the scan worker from the
components and RIs the search used, so it costs no second pass; a failed write, or a name two runs of the
batch map to (two `x.cdf` into one folder, neither is saved; names compared without case on Windows and
macOS), or a name that is itself an input of the batch, is an error of that run and keeps its result.
`.cmp` inputs are not saved again. In a folder a FastFlight run (`r1/MsData.FFT`) is saved as `r1.cmp`.
`write_cmp` writes a new `<name>.<pid>.part` and renames it, so a failed write leaves no broken `.cmp`; an existing
file is replaced, as by the old wizard.

Test `crates/md-core/tests/old_results.rs` checks files the old binary wrote itself (never committed):

```sh
# Linux/WSL with the old binary; writes example-data/old-results/ (ignored by git)
MD_RESULTS=$PWD/example-data/old-results docs/reverse-engineering/tools/golden/make-reference.sh
cargo test -p md-core --test old_results
# the old program reads what md-core wrote (the test saves it to old-results/md-core/)
docs/reverse-engineering/tools/golden/make-reference.sh --reread IN OUT [IN OUT]...
```

It compares every peak field with the old getters (9 digits), the spectra exactly, and the written
payload byte for byte with the old one; it passes without checking anything when the files are missing.

## Golden test

`crates/md-core/tests/golden_old_md.rs` compares peak detection and deconvolution with the old
MetaboliteDetector 3.3 (Linux) on the 15 example runs of `example-data/`; it passes without checking
anything when they are missing. The reference `crates/md-core/tests/data/old-md-3.3.tsv` holds the old peaks
(ion, first, apex and last scan, S/N) of 7235eg04, eg20 and eg32 and the components (apex scan, ions,
heights, integrals) of all 15 runs.

How the reference is made: `docs/reverse-engineering/tools/golden/make-reference.sh` (Linux or WSL, old
`.deb` installed) loads a small driver into the old binary with `LD_PRELOAD`. The driver replaces `main`
and calls the binary's own functions in the order of the GUI's "Detect components": `CdfImporter`,
`GCMSDiskScan`, `GCMSAnalyzer::detectAllPeaks(true, 10, 0.05)`, `deconvoluteChromatogram`. It uses the
static `AN_*` defaults compiled into the binary, which are also the values of a fresh config file. This
is the old code itself, not a re-implementation. The GUI export was not used because it needs clicks
per run. The script checks the binary's SHA-256, so
another build cannot silently produce a different reference. Neither the binary nor the example data are
committed.

Tolerances and their reasons:

| Compared | Tolerance | Reason |
|---|---|---|
| Peaks: ion, first, apex, last scan | exact | same algorithm and integer data |
| Peak S/N | 1e-5 relative | the reference prints 6 significant digits |
| Component ions and apex scan | exact | same fit, bit for bit (`quadfit`) |
| Ion height | ±1 | the old code rounds a few modelled heights the other way |
| Ion integral | 1e-5 relative | the reference prints 7 significant digits |

All 246 components match. Until issue 250, 5 differed: at precise-apex ties the old GSL fit rounds by
its own rounding errors. `quadfit` now rebuilds that fit from black-box outputs only, so no GSL code is
copied (GSL is GPL, this project MIT); the method is in `docs/reverse-engineering/gsl-clean-room.md`.

A second reference, `old-md-3.3-matched.tsv`, is made the same way with the old `MatchedPeakDetector(10, 3)`
(`MD_MATCHED="10 3"` for the driver; the *Calibration Setter* values, §3.3) and compared with
`Detector::MATCHED` under the same tolerances. Its peaks (12 590, 12 411 and 12 422 in eg04, eg20, eg32)
match exactly, S/N included (the old default -1), and all 2 807 old components match. The filtered traces
are integers with broad, smooth tops, so the parabola vertex lands exactly on a half bin much more often
(61 of the 12 590 peaks of eg04) and flat tops are common; before issue 250, 150 components differed
because of that. Matched peaks can also start at their apex: the old program then leaves the precise apex
at 0, which puts the peak into histogram bin 0 (10 components of eg07–eg25 depended on it).

A third reference, `old-md-3.3-high-resolution.tsv`, covers a synthetic high-resolution run made from 7235eg04:
- every m/z moves up by 0.047;
- a second fragment at 147.089 repeats the 147 trace 5 scans later.

To regenerate it, run `make-reference.sh --high-resolution`. The run itself (`example-data/high-resolution/`) is not committed. The test builds the run in memory from eg04, and the script writes it as a CDF for the old app with an ignored test (`write_high_resolution_run`).

The old CDF import centroids such scans to nominal masses (§4.1). The port does the same, so all peaks and components match with no known difference.

A fourth reference, `old-md-3.3-exact-masses.tsv`, takes the same run as a disk-scan cache with DataType 1, as FastFlight writes it. The old app then keeps the exact masses and runs `reorderPeaks`, so 147.047 and 147.089 stay apart in the components. The port matches it with no known difference. The cache carries per-scan statistics with nf = 1 and mean = sd = 0 (§2.3), not the `analyzeBaseline` values of a real FastFlight import. To regenerate it, run `make-reference.sh --exact-masses`; the ignored test `write_exact_mass_run` writes the cache.

The test fails when one of these counts grows.

### Stage fixtures

`crates/md-core/tests/fixtures/legacy/RUN.{peaks,regions,components}.tsv` hold the old program's
intermediate output for each of the 15 runs (default detector, static `AN_*` defaults): every peak of
`detectAllPeaks` in list order, every region of `detectComponents` (§4.5, rejected ones too, as peak
lines), and every compound with its model peak and per-ion fit. `make-reference.sh --stages` writes them
(the same driver, plus `stages.py`; format there), about 1.3 MB. `tests/legacy_stages.rs` compares
`find_run_peaks`, `deconv::component_regions` and `analyze_run` against them, in order:

| Compared | Tolerance | Reason |
|---|---|---|
| Peak ion, first, apex, last scan; peak count and order | exact | same algorithm |
| Interpolated apex | 5e-4 scan | md-core keeps it as f32 |
| Peak S/N, quality, sharpness | 1e-6 relative | md-core computes them in f32 |
| Component count, order, apex; model peak bounds and values | exact | integer sums |
| Regions (`deconv::component_regions`): count, order, peak lines | exact | same algorithm; bin bounds not compared |
| Ion height, model share `c2` | exact | `legacy_numeric::round_sig8` (#365) |
| Ion integral | 5e-9 relative | the dump's 9 digits |

The golden test found three porting errors in peak detection, now fixed:

- After an accepted peak the old loop steps back one scan, so a tailing peak starts right after the
  previous one.
- The apex search starts at the last rising scan, so a tie keeps the later scan.
- The baseline fit takes half of the distinct residuals (`std::map` keys), not half of all points, and
  the pivot tie goes to the right.

## Libraries

`md_core::library` reads the old binary `.lbr` libraries (gzip or raw, versions 2 to 5, layout in
[file-formats.md](reverse-engineering/file-formats.md) §1), imports and exports AMDIS MSL (§6), and imports
NIST text (`.sdf`, §6.3) and NIST MSP text (`.msp`, §6.4). NIST data is licensed and never committed; the
tests use the synthetic `tests/data/synthetic-nist.sdf`.
Library files come from users, so reading is a trust boundary: every length is checked against the
remaining bytes before anything is allocated, check values are verified, gzip output is capped at 1 GiB,
and errors are returned instead of panics. Versions 4 and 5 carry peaks and a model peak per compound
(§1.6, `LibraryCompound::peaks` and `model`), version 3 adds `formulas`, version 5 `mid_sizes`. Version 2
and 3 compounds with peak records are rejected as not supported (their layout is unknown). Strings are read
as UTF-8, falling back to Latin-1. The 40-byte peak records of versions 2 and 3 (§1.3) are kept verbatim
as `blocks`. `tests/lbr_sample.rs` checks the sample library `lib_aa.lbr` from the old package when
`MD_LBR_SAMPLE` points to it, and with `MD_LBR_V5_DIR` the version 5 files that the old program wrote and
read back (`tools/golden/make-lbr-v5.sh`); none of these files is committed.

### Library text format (.mdlib)

Edited libraries are saved in a new text format, not as `.lbr`: a text file can be read, diffed and
repaired by hand. It needs no new dependency (no
serde in `md-core`), and unlike JSON it keeps NaN and infinity. `write_mdlib` and `parse_mdlib` keep
every field of `LibraryCompound` bit for bit (only a NaN loses its sign and payload), so `.lbr` →
`.mdlib` → reading back loses nothing.

```
MetaboliteDetector library 2
compound
name Serine 3TMS
scan 0
rt_ms 436080
ri 1369.05
bin_width 1
peaks 73:1000 204.1433:250
fwhm null
integrals null
quant_ions 204.1433
auto_quant_ions
feature FRAG	(204/3)(218/2)
feature COMMENT	from AMDIS
drt 0
time_stamp 1227546843
id 2197790076
comb_spec 73:1000
ion_freq 73:1
spec_count 1
ris 1369
block 0707…07
source Initial import from/home/kh/Desktop/test.msl
peak 204 1 250 0 -1 1 2169 2363 1 2223 0
model 73 -1 -1 0 -1 0 0 0 1 0 0 0.5 0
formula 73 147 C3H9Si
mid_size 73 147 2
end
```

- UTF-8 text (a leading byte order mark is skipped), one `key value` per line (value = everything after the first space). Blank lines and lines
  starting with `#` are ignored. The first other line is the header; its number is the format version.
  A reader rejects a newer version and every unknown key, so a new field means a new version. Version 2
  added `peak`, `model`, `formula` and `mid_size`; version 1 files still read.
- `compound` … `end` is one compound. Inside, `name` is required; every other key is optional and
  defaults to a new compound's value (as after MSL import). Each key appears at most once, except
  `feature`, `block`, `source`, `peak`, `formula` and `mid_size`, which repeat in order.
- Numbers are written as Rust prints them: the shortest text that reads back to the same `f32`/`f64`
  (`NaN`, `inf`, `-0` included). Lists are space-separated; `m/z:value` pairs are joined by `:`.
  `peaks` holds the spectrum (`masses`, `values`); `fwhm` and `integrals` are `null` or a list (the old
  null pointer differs from an empty list). `quant_ions` is sorted and deduplicated on reading.
- Text (`name`, feature key and value, `source`) escapes `\\`, `\n`, `\r` and `\t`; any other escape is an
  error. A feature is `KEY<tab>value` with a unique key (MSL import makes an empty one from `: x`). MID fragments are the feature `FRAG`
  (syntax `(174/8)(261/5)`, F16) and comments the feature `COMMENT`; both keys occur in the old binary.
- `block` is one 40-byte record of `.lbr` §1.3 as 80 hex digits.
- `peak` is one `.lbr` §1.6 peak: ion, sharpness, integral, quality, signal to noise, count, first scan,
  last scan, model parameter, max scan, interpolated max scan. `model` is the same 11 fields, then the
  noise factor and the model values. `formula` is `begin end formula` (formula escaped like text, may be
  empty), `mid_size` is `begin end size`.
- Reading is a trust boundary like `.lbr`: files over 1 GiB are rejected, every error names its line.
  `save_mdlib` writes `<file>.tmp`, syncs it to disk and renames it over the target, so a failed save
  never leaves half a library.

The library editor (File → Library editor…, `crates/md-app/src/library_editor.rs`) edits a copy of
the open library or a new one: create, delete and rename compounds, add the selected component of the
active run ("Add component…", see below), copy RT and RI of the run's assigned components ("Copy RTs…"), view the spectrum, and edit RT, RI, quant ions, MID fragments,
comments and other features. Collapsed sections (`library_editor/extras.rs`) show the automatic quant ions
("Use as quant ions", "Clear automatic ions"), the ion frequency with the spectrum count (read only), the resolution
(`bin_width`, read only) and the v5 MID fragments as a table (start, end, sum formula checked with `isotope::parse_formula`, MID size:
add, change, remove). That table is the list the MID wizard reads (`mid::series::fragments`); the `FRAG` feature is a separate list
that only the FRAG export reads. It saves only `.mdlib`. Save writes back only a `.mdlib` file it opened or
saved before; every other target goes through the native save dialog, which asks before replacing a file,
and a target that is not `.mdlib` is refused, so `.lbr` and `.msl` files are never overwritten. "Export MSL" writes the library as AMDIS text in the three
types of the old export dialog (Nominal, Profile, FWHM; `library::write_msl_as`, file-formats.md §6.2) and only to a
`.msl` file. After a
save the app opens the saved library, so the next search uses the edits.

**Growing a library from runs** (`md_core::library_add`, old `DialogAddLib` and `copyRT`, algorithms.md §7e).
"Add component…" turns the selected component into a compound (deconvoluted spectrum, RT at the apex, RI
from the run's calibration or -1, automatic quant ions by §7.2) and puts it into the edited library as a
new entry, or, for a name already in the library and after a yes/no question, replaces the entry, replaces
only spectrum, RT and RI (quant ions and comments stay), or sets only the automatic quant ions. The quant
ions come from the peaks the run's components were built from (`Run::peaks`, also those of a loaded
`.cmp`), so they do not depend on the current deconvolution width. "Copy RTs…" searches the components of the active
run against the edited library, lists the assigned ones with their old and new RT and RI and, on "Copy",
writes RT and RI into those entries. Both change only the copy in the editor until Save; "Modify entry"
(spectrum merge) is not ported.

## Old raw-data cache (.bin + .idx)

The old app converted every imported NetCDF into `<run>.bin` (scan records) and `<run>.idx` (index);
some users only kept these. `md_core::diskscan::read` reads the pair into the same `Chromatogram` as
`andi::read` (layout in [file-formats.md](reverse-engineering/file-formats.md) §2). Reading is a trust
boundary like the libraries: every offset and count is checked against the file size before anything
is allocated, each `.bin` record header must match its index entry, and errors are returned instead of
panics. Caches with a TOF calibration polynomial are rejected (the old app applied it to the masses;
`CdfImporter` writes none).

`andi::write` exports a `Chromatogram` as ANDI-MS NetCDF (64-bit-offset NetCDF-3): the raw-data
variables `scan_acquisition_time`, `scan_index`, `point_count`, `mass_values`, `intensity_values` and
`total_intensity` (the sum of each spectrum). `md_core::read_run` picks the reader by extension, so the
app's File menu and the batch view open `.bin` files; File → Export as NetCDF writes the active run.

Proof: `docs/reverse-engineering/tools/diskscan/make-cache.sh` converts the 15 example runs with the
old binary's own `CdfImporter` (`LD_PRELOAD` driver, like the golden test). The test
`old_cache_matches_cdf` in `tests/example_data.rs` (runs when `MD_OLD_CACHE` points to that output)
checks that every cache has the same spectra as the original `.cdf`, and that its NetCDF export reads
back unchanged. Times match to f32 precision: the old importer kept times as f32 seconds. The old binary
also re-imports the exported NetCDF files with identical spectra and times; only the stored TIC differs,
because the example `.cdf` files carry the instrument's `total_intensity`, not the sum of the spectrum.

## Spectrum normalisation and export

`md_core::spectrum::normalize` divides a spectrum by one norm: none, base peak, unit (vector length) or
relative (total signal), like the old `SpectrumPlot` ([algorithms.md](reverse-engineering/algorithms.md) §7d).
Spectra are centroid in this rewrite (NetCDF and `.bin` hold ions, not profile points), so there is no
profile/centroid switch. File → Export (`crates/md-app/src/export.rs`) writes the shown chromatogram of the
active run (TIC and/or SIC, every scan, not the decimated plot points) and the shown spectrum (with the
chosen normalisation) as CSV through `md_core::csv::write_columns`, and each plot as PNG. PNG is a crop of
egui's own window screenshot, so it looks exactly like the screen. Files are written to `<name>.part` and
renamed.

## Component panel and spectrum comparison

`md_core::detail` holds the values of the old component panel (mean S/N, mean sharpness, peak discrepancy
index, the peak list with one row per ion) and `compare`, which puts the found and the library spectrum on
one mass axis, normalised with `spectrum::normalize` ([algorithms.md](reverse-engineering/algorithms.md)
§7g). `crates/md-app/src/detail.rs` draws the right panel of the selected component and the mode picker
above the spectrum plot: found only, split (library mirrored below) and difference (found minus library).
A click on an ion in the peak list or on the spectrum sets the SIC to that mass. The compared library hit is
the assigned one, else the best. "Copy both" and "Save both as CSV…" export both spectra and their difference.
`Run` keeps the `Peak` list of the analysis so the panel can read S/N, quality and sharpness.

## Biomarker statistics and PCA

`md_core::stats` works on the batch table as shown (after normalisation, so the same numbers as the
Batch window). `compare` takes two lists of replicate groups and gives, per substance, hits, mean, SD,
ratio and a two-sided p-value: Welch's t-test over all runs (old "independent"), over one mean per
replicate group (old "not independent"), or a paired t-test (new). A missing value counts as 0 by
default, like the old program. The t distribution is our own regularised incomplete beta function; no
statistics crate. `five_numbers` gives min, quartiles (R type 7) and max for the box plots. `pca`
excludes the pool runs, z-scores each substance, and takes the SVD of the runs × substances matrix with
`nalgebra` (pure Rust, `std` only, no BLAS). This equals the old covariance eigen-decomposition and
scales with the number of runs squared, not the number of substances squared, which matters for
non-targeted tables. The sign of a component is chosen so its largest loading is positive. Values are
tested against R 4.5.2 (`scripts/r/stats_reference.R`, relative 1e-9). Details:
[algorithms.md](reverse-engineering/algorithms.md) §7f. The window
(`crates/md-app/src/statistics.rs`) recomputes only when the batch table or the groups change; which rows it lists and how numbers read comes from `md_core::view`.

## Mass isotopomer distributions (MID)

`md_core::mid` follows the old MID wizard (algorithms.md §7g). `mid_library` keeps the library substances
with fragments (sum formula per ion range) and makes every fragment ion a quant ion; the app scans the
chosen runs against it like a batch, with the rescan filling fragment ions that were not deconvoluted,
and `determine` turns the integrals into MIDs. Uncorrected and reference correction give one MID per
fragment over the labelled (and unlabelled) runs; library spectrum and sum formula give one per run.
Substances are solved in parallel (rayon); `examples/perf.rs` times 1 000 substances × 200 runs. The batch
window's "Export MIDs (FRAG)" writes the old export of the `FRAG` feature, e.g. `(174/8)(261/5)`, from the
shown table (`mid::frag_export`). The worker also computes the mean and standard error of the labelled
runs (`mid::group_stat`, per-run methods) and the old result window's text exports (`Mids::run_text`,
`stats_text`, `group_text`, numbers in Qt's `%g` via `mid::qt_g`); the window shows them as table
columns and bar plots and saves or copies the texts. ‹ › move a run column and ▲ ▼ a substance with its
tables and plots; the window then swaps the fits or substances in its `Mids` and rebuilds the statistics
and texts, so the exports follow the shown order (old reorder dialogs and deviations: algorithms.md §7g.7). Methods 2 and 3, the statistics and both per-run
exports are tested against the old `MIDCompoundSeries`, methods 0 and 1 and the group export against the
old `MIDWizard::Calculate` and `MIDResultWidget` (`tests/mid_series.rs`, harness
`docs/reverse-engineering/tools/isotopes/make-midseries.sh`; the old program computes but never shows the
group results, algorithms.md §7g.2), the solver
against `tools::MIDSolver`, and the pipeline on the unlabelled example runs (M+0 ≈ 1).

## Settings

`md_core::settings::GcmsSettings` holds all old settings with the old defaults (file-formats.md §5.1).
`GcmsSettings::analysis()` derives `AnalysisSettings { peaks, deconv, quant, search, fwhm }` once; it is the only
conversion, and `ScanSettings` carries it. One const table, `settings::KEYS` (old key, display text (English, German),
accessor, range, fixed slot in the old file; no macro), drives the app's `KEY=value` lines (`set`, `lines`), the
old file's fixed-size fields, `validate` (`sanitized` resets the bad ones) and the ranges of the Settings window
(`limits`). Its keys are the constants of `settings::keys` (the app names no key as a string literal); a new value is
one constant, one table row plus its struct field and default. `an_number_bins` is read once in `analysis()`
for the peak and the deconvolution step. It reads and writes the old settings
file (§5.2: raw, 191 + 8n bytes, no header; a gzip wrapper is accepted because algorithms.md §1 reports
one). The file is untrusted: size must match the excluded-ion count exactly, the checksum must match,
every value must be finite and in range, and gzip output is capped at 1 MiB. A text variant is not read
(not decoded). `LS_CUTOFF`, `AN_MIN_HEIGHT_BP`, `AN_MODELED_VALUES_FOR_INTEGRATION` and `LS_MASS_FILTER`
are not in the old file; loading keeps their current value. `SF_CENTROID_*` steer the centroid export.
`SF_FWHM` (0..1) is the share of a peak's maximum a point must reach to count in every centroid: the
high-resolution import (`IonTraces::new`), the centroid export, the MSL export and the TOF calibration (#293).
`SF_PEAK_*`, `SF_MAX_BASE_DEV` and `LS_USE_SCALED_LIB` are dead in 3.3 too (algorithms.md §4b.4); they and
`AN_MAX_BASE_DEV`, `LS_RESOLUTION` and `QUANT_ION_NUMBER`/`_MIN_DIST_IONS`/`_QUALITY_CUTOFF` are
stored and shown, but no step reads them yet. The search cuts off at `LS_CUTOFF` everywhere, unlike the old GUI (§5.8, #399);
the component list, Find substance and the batch use the same cutoff. The Settings window also chooses the two libraries and holds the old Appearance switches (axis
titles, m/z labels on the 5 highest peaks in view), which apply at once.

In the app (`crates/md-app/src/settings.rs`, Analysis → Settings…) edits go to a draft until Apply.
Applying re-runs `IonTraces` → peaks → deconvolution for every open run on worker threads (stale jobs are
dropped), clears the selected component and library hits, and re-prepares the library when `LS_*` changed.
New runs open with the applied settings. `Settings::analysis` is the applied `analysis()`, derived once when the settings are applied (`set_applied`); the UI reads it per row and frame instead of deriving it again (#415).

## Saved state

`crates/md-app/src/persist.rs` (S4, #253). eframe's `persistence` feature writes one file, `app.ron` in the
config dir (below), on exit and every 30 s: egui's memory (panel sizes, tool window positions), the window
geometry and, under the key `metabolite-detector`, our `KEY=value` lines. These hold every `GcmsSettings`
field under its old key (`LS_RI_DIFF=10`, ranges as `73:75;147:148`, `PEAK_DETECTOR=default` or
`matched <width> <threshold>`), `panel.<area>` for the View menu, `library`, `nist_library` and `recent`
(newest first, at most 10; File → Open recent) and `mode` (`guided` or `expert`, S7). Keys the version does
not know are written back unchanged. `language` is `en` or `de` (S8), saved only once the user picks a
language other than the system's; without it the app follows the system language.

## UI language

S8, #257. `md-app/src/lang.rs`: every UI text is an English/German pair at its call site, `lang::t` for
fixed texts and `lang::tf!` for formatted ones, so no text can miss a language. The language is per
thread (a thread starts in English, so parallel tests do not see each other's choice); `main` sets the
system language (`sys-locale`), `persist` the saved one, and `Job::spawn` hands it to its worker. CSV
headers, file contents, default file names and saved keys stay English literals; where md-core's English
text is shown (`Normalization::label`, "No match: N") md-app translates at the display site. A window
whose title is translated keeps the English title as its id. `every_ui_text_is_translated` scans the
sources for bare widget texts.

md-core errors (#339): md-core, md-cli and files keep the English `Display` texts. `md-app/src/errors.rs`
translates them where the window shows them: `errors::text` matches a message against a table of
English templates (`errors.txt`, `{}` for the changing parts, translated in turn; a file name before `: ` stays) and
returns the English text when none fits (an operating system text). A red error label draws a fixed text,
a `Msg` or `errors::text` (`error_labels_translate_core_messages`), and `real_core_messages_are_translated` runs real md-core errors
through the table, so a changed wording of one of those texts fails there. A new md-core error text needs a template.

Stored texts (#341): a note, an error or the RI source that is shown later is a `lang::Msg`, not a `String`.
`Msg::new(|| ...)` (or `msg!`, `tf!` as a `Msg`) runs its text once per language and keeps both; the language is
picked when it is drawn (`&*msg`), so a switch changes texts already on screen. `md-core` texts and file names
enter through `Msg::from(String)`, the German side through `errors::text`. Jobs build a `Msg` too, so the language
of the worker thread does not matter (`jobs::poll_job` returns `Result<T, Msg>`). egui_plot hides a curve by its
id, which is its name unless set: a curve with a translated name needs a fixed `.id(...)` when its plot has a legend.

## Guided mode

`crates/md-app/src/guided.rs` (S7, #256). Without a saved `mode` the start asks (a modal, Guided
preselected); the answer is saved, and Help and Settings switch it. The mode is published to the egui context
each frame; `guided::advanced` puts expert options of the Settings, RI calibration, batch and MID windows in
an "Advanced" section, closed in Guided and open in Expert mode (salted with the mode, so a switch resets it).
Analysis → Step by step starts one of four flows (RI calibration, component detection, batch quantification,
MID) in its own window: "Step x of n", Back, Next (once the step's input is there), Skip, and Finish or "Open
the window". Each step draws a building block of the expert window (`Calibration::contents`, `Batch::runs_row`,
`MidWindow::settings`, …) on the same state, so the flow ends in that window with its values set. Every
string of the module is in `guided::text` for the translation (S8).

Every value is checked at start: a bad value keeps its default, and each field `GcmsSettings::validate`
names falls back (`sanitized`, in md-core) to its default until the settings pass. Then, or when eframe could not parse `app.ron` (our
key is missing from an existing file), the file is copied to `app.ron.<unix s>.bak` and the Settings window
opens with the message. A library that is missing at start (a drive not mounted) is not opened but kept.
Paths are saved absolute; a recent run that is clicked leaves the list until it loads again.

First start (no `app.ron`): the old 3.3 settings are imported once and the Settings window says so. Sources:
`HKCU\Software\TU-Braunschweig\MetaboliteDetector` (through `reg query`) on Windows, else
`$XDG_CONFIG_HOME` or `~/.config` + `/TU-Braunschweig/MetaboliteDetector.conf` (Qt INI: quoting,
`\xHH`/`\0` escapes; `LS_RESOLUTION` is a `@Variant` float, type 135, big-endian; Qt on Windows writes ints as
`REG_DWORD` and that value as `REG_BINARY` UTF-16). A saved value wins over the default (`LS_RI_DIFF` 10 over
20). `AN_BP_HEIGHT` and `MASS_FILTER` are read as `AN_MIN_HEIGHT_BP` and `LS_MASS_FILTER` (Inferred: same
meaning and default); non-empty old range lists are not decoded. `rc_lib_file`, `rc_lib_nist_file` and
`rc_open_file` become the libraries and the first recent run, `[GUI] size` the window size; other
`[Algorithm]` keys (`baseline_correction`, `rc_*`) are kept for later versions, `LS_CUTOFF` and
`[GUI] axes_label`/`peak_label` are read, the `[GUI]` window state blobs are not.

## Updates

`crates/md-app/src/update.rs` (Help menu, #235): `update/release.rs` holds the pure parts (newest
non-draft release above the running version from the GitHub releases JSON, `SHA256SUMS.txt` lookup,
asset name per install kind, install kind from the exe path, `update.txt` settings), `update/install.rs`
the network and file work. Release builds check at start when enabled and the last check is a day old;
debug builds never check on their own and never install. Settings: `update.txt` (`auto=1`,
`last_check=<unix s>`) in `%APPDATA%\MetaboliteDetector` or `$XDG_CONFIG_HOME/metabolite-detector` (else `~/.config/metabolite-detector`).

After the click the app downloads `SHA256SUMS.txt` and the one asset named for its platform and version
over HTTPS, and installs nothing unless the SHA-256 matches:

- MSI (exe in `%LocalAppData%\Programs\MetaboliteDetector` or, with `ALLUSERS=1`, `%ProgramFiles%\MetaboliteDetector`):
  saved to a fresh temp dir; a hidden PowerShell (`-EncodedCommand`, paths as environment variables) waits for the
  app to exit, runs `msiexec /i … /passive /norestart`, removes the temp dir and starts the exe again.
- Portable exe (anywhere else on Windows): written as `<exe>.new`, the running exe renamed to `<exe>.old`, the new
  one moved in and started; the next start deletes `<exe>.old`.
- Linux tar.gz: `metabolite-detector` (the running exe) and `md-cli` next to it, if present, are written as
  `<file>.new` with the mode of the file it replaces and renamed over the originals only once all are written.

## Build and test

```sh
cargo build --release            # binary: target/release/metabolite-detector(.exe)
cargo test --workspace
cargo clippy --workspace --all-targets
cargo run --release -p md-core --example perf
cargo run --release -p md-core --example batch -- 4 library.msl run1.cdf run2.cdf
cargo run --release -p md-cli -- --library library.msl -o table.csv run1.cdf run2.cdf
cargo run --release -p metabolite-detector -- path/to/run.cdf
example-data/fetch.sh             # Windows: .\example-data\fetch.ps1
```

In a git worktree, `node scripts/link-example-data.mjs` links the main checkout's example data (`.ps1` and `.sh`
do the same), so the golden tests run instead of skipping.

`.github/workflows/ci.yml` runs for every non-draft PR and every push to `main`, with no paths filter, because the required
checks must report on every PR. Ubuntu runs `cargo fmt --check`, Clippy (warnings fail), the tests through `cargo llvm-cov` and the
Sonar scan, so Linux builds once for both; Windows only runs the tests, since no code is platform-specific. PRs without changes to
`crates/`, the Cargo files, the toolchain, `sonar-project.properties` or the workflow skip the build and only run the Sonar scan;
`main` always builds. CI builds without debug info. The summary job `CI` fails unless both OS jobs succeeded; it is the required CI
result for `main`. It does not fetch the example data, so those tests skip themselves;
run them locally. Only `main` saves the cargo cache; PRs restore it. The cache key
follows `rust-toolchain.toml` and `Cargo.lock`, not the toolchains a runner image happens to preinstall
(CI uninstalls those before restoring the cache).

## Compatibility goals

- Read the old binary files students already have: `.lbr` libraries first, then `.cmp`, `.calib`,
  `.bin`/`.idx` (see [reverse-engineering/file-formats.md](reverse-engineering/file-formats.md)).
- Reproduce the old numerical results (peaks, compounds, scores, areas, MIDs) within a documented
  tolerance, verified against output of the old app on the same input.
