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
| NetCDF | `netcdf3` (pure Rust) | ANDI-MS files are NetCDF-3; no libnetcdf/HDF5 C dependency to ship on Windows. |
| Parallelism | rayon | Data-parallel loops over scans, ions and files. |
| R console | Not ported directly | RInside tied the old app to a Linux R install. Replace with CSV/R-script export first; optional embedded scripting later. |

## Layout

```
Cargo.toml              workspace, release profile (LTO, 1 codegen unit)
src/md-core/            library: data model, import/export, analysis algorithms (no GUI deps)
  src/lib.rs            Chromatogram (flat ANDI-style arrays), TIC/SIC, plot decimation
  src/andi.rs           ANDI-MS NetCDF import with input validation
  src/library.rs        compound libraries: old .lbr (version 2) reader, AMDIS MSL import/export
  src/peaks.rs          ion traces (ion-major index), noise, baseline and peak detection per trace
  src/deconv.rs         deconvolution of ion-trace peaks into components with spectra
  src/ri.rs             retention index: RT ↔ RI calibration, marker search, RT shift, .calib file
  examples/perf.rs      performance smoke test on synthetic large runs
src/md-app/             egui desktop app (binary `metabolite-detector`), drawing only
  src/calibration.rs    RI calibration window (alkane search, .calib load/save)
docs/reverse-engineering/  evidence from the old binary (features, algorithms, file formats)
```

Planned additions: `src/md-cli` (headless batch runs for big studies) and modules in `md-core` for
library search, quantification, MID and statistics.
All analysis code lives in `md-core` so the GUI, CLI and tests share it.

## Performance rules

Large studies (hundreds of runs, tens of millions of points each) are the normal case.

1. **Compact storage.** Spectra are flat `f32` arrays with scan offsets (≈ 8 bytes per data point),
   matching the old `GCMSDiskScan<int, float>`. Accumulate sums in `f64`.
2. **Never block the UI.** File loading and every analysis step run on worker threads; the UI polls
   results and shows progress. Long jobs must be cancellable.
3. **Draw proportional to pixels, not data.** Line plots go through `md_core::decimate`
   (min/max per pixel bucket of the visible range), so peaks stay visible at any zoom.
4. **Parallel by default.** Per-scan, per-ion and per-file work uses rayon.
5. **Index once, query many.** m/z is sorted within each scan (binary search for SICs). Peak detection
   and deconvolution share one ion-major index per run (`IonTraces`) instead of rescanning spectra per ion.
6. **Measure.** `examples/perf.rs` reports timings; each new algorithm adds a case, and PRs that
   touch hot paths report before/after numbers against the table below.
7. Memory-mapped or streamed access (like the old disk scan) is added when real runs exceed RAM;
   the cut-off is recorded as a `shortcut:` comment in `andi.rs`.

Measured 2026-10-19, release build, 36 000 scans × 500 ions (18 M points, 137 MiB):

| Operation | Time |
|-----------|------|
| TIC | 2.0 ms |
| SIC (one m/z) | 1.0 ms |
| Decimate TIC for drawing | 0.05 ms |
| Overlay frame: decimate TIC + SIC of 15 runs | 0.8 ms |
| Ion traces (index, 500 traces) | 26 ms |
| Peaks in all 500 traces (300 000 peaks) | 20 ms |
| Deconvolution (12 000 substances of 25 ions) | 19 ms |

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
detector on every trace in parallel. `PeakSettings::default()` holds the old `AN_*` defaults.
`detect_peaks_matched` and `matched_filter` port the old `MatchedPeakDetector` and its filter, which the
old app only uses in the RI calibration dialog (and deconvolution, §4.5).

![Detected peaks on four traces of 7235eg04.cdf](images/peaks-7235eg04.svg)

Conscious deviations from the old app, each to be checked against old output in the golden-test issue:

- **Traces per nominal mass.** m/z is rounded and points of one nominal mass in a scan are summed. The
  old app's trace construction for ANDI data is unknown; the example runs are unit-resolution, so this
  matches them. High-resolution data needs exact-mass traces and `reorderPeaks` (§4.2).
- **Noise factor 1.0.** The old detector keeps its constructor value 1.0; `GCMSAnalyzer` overwrites it
  from a source not found yet (§2.3). The noise factor and MAD of `analyzeBaseline` are not computed
  because the detector never reads them.
- **f32 intensities instead of int.** Only the SG smoothing and the matched filter output are truncated
  to integers like the old int filters. S/N is not rounded to 8 significant digits.
- **No progress or cancel per 100 traces.** A whole run takes milliseconds (table above).
- **Border extension off** like the old default (`setPeakBorderExtensionEnabled` has no known caller),
  so the extension code is not ported.

## Deconvolution

`md_core::deconv::deconvolute` ports the old `deconvoluteChromatogram` ([algorithms.md](reverse-engineering/algorithms.md)
§4.4–4.9): the apex intensities of all peaks go into a histogram of `AN_NUMBER_BINS` bins per scan, the
matched filter (`AN_DECONVOLUTION_WIDTH` scans) marks regions above 5, and each region becomes a component.
Its model peak is the sum of the sharpest clean peaks (quality ≤ 10, at most a quarter); every ion trace is
fitted as `c0 + c1 * scan + c2 * model`, and `c2 * max(model)` is that ion's height in the spectrum.
Components are fitted in parallel. `DeconvSettings::default()` holds the old `AN_*` defaults. The app runs
peak detection and deconvolution after loading and marks each component on the TIC; a click near a dot
shows its deconvoluted spectrum.

The old quirks are kept on purpose until the golden test (issue 17) has compared both versions: the first
and last bin of a region are added twice, so a peak can appear twice in a component and counts twice
towards `AN_MIN_PEAK_NUMBER`; the model peak sums raw intensities including the baseline, so all heights
of a spectrum share one scale factor.

Conscious deviations:

- **No precise apex → apex scan.** Where the parabola fit was refused, the histogram uses the apex scan.
  What the old field holds then is unknown.
- **Collinear fits drop the peak.** The old code uses an SVD (`gsl_multifit_linear`), which also returns a
  minimum-norm fit when model and scan are collinear; here those peaks are skipped.
- **Filter edges are zero-padded** (old edge handling unknown, §4.11), and `roundSignificantFigures` is
  re-implemented (round to 8 significant digits).
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
(m/z 57, C10–C36 by default) or loads/saves a `.calib`; the calibration applies to all open runs, and the
status bar shows the RI of the scan or selected component.

Checked against TargetSearch on the example data (`tests/example_data.rs`): its FAME markers (m/z 87 in
the windows of `rimLimits.txt`) give the RI of every scan in `RI_*.txt` to within 0.01 in 14 of 15 runs.
In `7235eg06` TargetSearch took a peak at 321.1 s as the middle marker, where m/z 87 reaches 36 512;
the marker peak at 311.1 s, as in all other runs, reaches 221 839. Filter width 40 scans, threshold 100.

RT shift correction (§6.5): `pick_reference` chooses the substance with the best library score for the
reference near its RI, `shifted_to` shifts the calibration so it lands on the library RI. The scores
come from library search (issue 19); the app does not offer the correction yet.

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

## Libraries

`md_core::library` reads the old binary `.lbr` libraries (gzip or raw, version 2, layout in
[file-formats.md](reverse-engineering/file-formats.md) §1) and imports and exports AMDIS MSL (§6).
Library files come from users, so reading is a trust boundary: every length is checked against the
remaining bytes before anything is allocated, check values are verified, gzip output is capped at 1 GiB,
and errors are returned instead of panics. Version 5 (the old writer's current version) and compounds
with peak records are rejected as not supported until a sample file shows their layout. Strings are read
as UTF-8, falling back to Latin-1. The unknown 40-byte block list (§1.3) is kept verbatim.
`tests/lbr_sample.rs` checks the sample library `lib_aa.lbr` from the old package when
`MD_LBR_SAMPLE` points to it; the file is never committed.

## Build and test

```sh
cargo build --release            # binary: target/release/metabolite-detector(.exe)
cargo test --workspace
cargo clippy --workspace --all-targets
cargo run --release -p md-core --example perf
cargo run --release -p metabolite-detector -- path/to/run.cdf
example-data/fetch.sh             # Windows: .\example-data\fetch.ps1
```

## Compatibility goals

- Read the old binary files students already have: `.lbr` libraries first, then `.cmp`, `.calib`,
  `.bin`/`.idx` (see [reverse-engineering/file-formats.md](reverse-engineering/file-formats.md)).
- Reproduce the old numerical results (peaks, compounds, scores, areas, MIDs) within a documented
  tolerance, verified against output of the old app on the same input.
