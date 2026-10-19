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
  src/andi.rs           ANDI-MS NetCDF import with input validation, NetCDF export
  src/diskscan.rs       old raw-data cache (.bin + .idx) reader
  src/library.rs        compound libraries: old .lbr (version 2) reader, AMDIS MSL import/export
  src/peaks.rs          ion traces (ion-major index), noise, baseline and peak detection per trace
  src/deconv.rs         deconvolution of ion-trace peaks into components with spectra
  src/ri.rs             retention index: RT ↔ RI calibration, marker search, RT shift, .calib file
  src/search.rs         library search: spectrum and RI scores, one-to-one assignment of components
  examples/perf.rs      performance smoke test on synthetic large runs
src/md-app/             egui desktop app (binary `metabolite-detector`), drawing only
  src/calibration.rs    RI calibration window (alkane search, .calib load/save)
  src/identify.rs       library loading, hit labels on the TIC and the hit list of a component
docs/reverse-engineering/  evidence from the old binary (features, algorithms, file formats)
```

Planned additions: `src/md-cli` (headless batch runs for big studies) and modules in `md-core` for
quantification, MID and statistics.
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
| Library search, 1 000 components × 10 000 entries, cutoff 0.8 | 31 ms |
| Library search, same, cutoff 0 (9 M hits kept, like the old GUI) | 250 ms |

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

Conscious deviations from the old app; the golden test (below) shows that none of them changes the
peaks of the example runs:

- **Traces per nominal mass.** m/z is rounded and points of one nominal mass in a scan are summed. The
  old app's trace construction for ANDI data is unknown; the example runs are unit-resolution, so this
  matches them. High-resolution data needs exact-mass traces and `reorderPeaks` (§4.2).
- **Noise factor 1.0.** `GCMSAnalyzer` sets the detector's noise factor from
  `DataProvider::getGlobalNoiseFactor` (§2.3). On all 15 example runs the old app logs
  "Error determining nf. Using nf=1", so 1.0 is what it uses there. Data where the old app finds another
  noise factor would differ. The noise factor and MAD of `analyzeBaseline` are not computed because the
  detector never reads them.
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

The old quirks are kept on purpose; with them, 225 of 246 old components match exactly (golden test below): the first
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
come from `md_core::search`; the app does not offer the correction yet (issue 45).

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
components in parallel and then assigns library compounds greedily by best overall score, each library compound
and each component at most once. `SearchSettings::default()` holds the old `LS_*` defaults.
With nominal masses each pair is one pass over the library peaks with a dense lookup into the component
spectrum, so no fragment pre-filter is needed (table above).

The app opens a library with File → Open library (or a `.lbr`/`.msl` path on the command line), identifies every
open run on worker threads, writes the assigned name at each component on the TIC (green above 0.9, orange
above 0.7, red otherwise, like the old app) and lists the best hits of the selected component.

Conscious deviations:

- **No RI yet → spectrum only.** The app does not pass the RI calibration to the search yet (issue 48). When the component or
  the library entry has no RI, the overall score is the spectrum score alone (with `use_spec` off there is then
  no hit). The old code then compares retention times / 1000 with the same sigma instead; that is not ported,
  because times of different runs and instruments only become comparable after calibration.
- **Nominal mass.** Library spectra are rounded to nominal mass and intensities of one mass are summed, because
  the old `getBinnedCentroidCompound` is not decoded (§5.9) and components are unit-mass. Matching is then exact,
  which makes the old 0.01 m/z tolerance and its double-match quirk (§5.3) irrelevant. Negative intensities
  and masses below 0.5 are dropped.
- **Hits need a score above zero.** The old search keeps pairs with score 0 and may even assign them; here they are
  no hits. The cutoff is applied during the search; the app searches with 0 like the old GUI. A higher cutoff
  yields the same assignments above it, because the greedy assignment takes pairs in descending score.
- **Ties are deterministic**: lower library index, then lower component index (old: pointer order).
- `LS_USE_SCALED_LIB`, `LS_IDENTIFICATION_ALGORITHM`, the Wave score variants and the SQL fragment pre-filter are
  not ported (semantics unknown, §5.9).

## Golden test

`src/md-core/tests/golden_old_md.rs` compares peak detection and deconvolution with the old
MetaboliteDetector 3.3 (Linux) on the 15 example runs of `example-data/`; it passes without checking
anything when they are missing. The reference `src/md-core/tests/data/old-md-3.3.tsv` holds the old peaks
(ion, first, apex and last scan, S/N) of 7235eg04, eg20 and eg32 and the components (apex scan, ions,
heights, integrals) of all 15 runs.

How the reference is made: `docs/reverse-engineering/tools/golden/make-reference.sh` (Linux or WSL, old
`.deb` installed) loads a small driver into the old binary with `LD_PRELOAD`. The driver replaces `main`
and calls the binary's own functions in the order of the GUI's "Detect components": `CdfImporter`,
`GCMSDiskScan`, `GCMSAnalyzer::detectAllPeaks(true, 10, 0.05)`, `deconvoluteChromatogram`. It uses the
static `AN_*` defaults compiled into the binary, which are also the values of a fresh config file. This
is the old code itself, not a re-implementation. The GUI export and `.cmp` files were not used because
they need clicks per run and `.cmp` holds no peak bounds. The script checks the binary's SHA-256, so
another build cannot silently produce a different reference. Neither the binary nor the example data are
committed.

Tolerances and their reasons:

| Compared | Tolerance | Reason |
|---|---|---|
| Peaks: ion, first, apex, last scan | exact | same algorithm and integer data |
| Peak S/N | 1e-5 relative | the reference prints 6 significant digits |
| Component ions and apex scan | exact for 231 of 246 | 15 known differences, issue 46 |
| Ion height | ±1 | the old code rounds a few modelled heights the other way |
| Ion integral | 1e-5 relative | the reference prints 7 significant digits |

Remaining differences, all in deconvolution (issue 46):

- 3 old components are missing.
- 12 components differ by 1–4 ions.
- 6 components have the same ions, but all heights and integrals differ by about 0.1–5 %. A different
  model peak is the likely cause.

The test fails when one of these counts grows. The counts are lowered as issue 46 fixes them.

The golden test found three porting errors in peak detection, now fixed:

- After an accepted peak the old loop steps back one scan, so a tailing peak starts right after the
  previous one.
- The apex search starts at the last rising scan, so a tie keeps the later scan.
- The baseline fit takes half of the distinct residuals (`std::map` keys), not half of all points, and
  the pivot tie goes to the right.

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
`total_intensity` (the sum of each spectrum). The app's File menu opens `.bin` files and exports the
active run.

Proof: `docs/reverse-engineering/tools/diskscan/make-cache.sh` converts the 15 example runs with the
old binary's own `CdfImporter` (`LD_PRELOAD` driver, like the golden test). The test
`old_cache_matches_cdf` in `tests/example_data.rs` (runs when `MD_OLD_CACHE` points to that output)
checks that every cache has the same spectra as the original `.cdf`, and that its NetCDF export reads
back unchanged. Times match to f32 precision: the old importer kept times as f32 seconds. The old binary
also re-imports the exported NetCDF files with identical spectra and times; only the stored TIC differs,
because the example `.cdf` files carry the instrument's `total_intensity`, not the sum of the spectrum.

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
