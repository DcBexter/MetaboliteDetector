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
  examples/perf.rs      performance smoke test on synthetic large runs
src/md-app/             egui desktop app (binary `metabolite-detector`), drawing only
docs/reverse-engineering/  evidence from the old binary (features, algorithms, file formats)
```

Planned additions: `src/md-cli` (headless batch runs for big studies) and modules in `md-core` for
deconvolution, library search, RI calibration, quantification, MID and statistics.
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
5. **Index once, query many.** m/z is sorted within each scan (binary search for SICs). Deconvolution
   will build an ion-major index once per run instead of rescanning spectra per ion.
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
| Ion traces (index, 500 traces) | 26 ms |
| Peaks in all 500 traces (300 000 peaks) | 20 ms |

Real runs (TargetSearchData, 4 400 scans, ~745 000 points, 6 MB `.cdf`): read + TIC 13–25 ms per
file; all 15 runs in 0.27 s; ion traces + peak detection 4–8 ms per file, 630–950 peaks in 236 traces
(`cargo test --release -p md-core --test example_data -- --nocapture`).
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
