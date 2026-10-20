# MetaboliteDetector 3.3 – Feature inventory

Source: the shipped Linux build `MetaboliteDetector-3.3.20200313-Linux_0.deb` (local `Old/`, not committed).
The ELF `usr/bin/metabolitedetector_app` (SHA-256 `7883d70c…4794`) is a Qt5/C++ GCC build with DWARF
debug info and full symbols. The source code is lost.

Evidence types: **S** = symbol/class names (`nm -C`), **T** = UI/log strings in the binary, **D** = Ghidra
decompilation, **P** = package contents. Algorithm details: [algorithms.md](algorithms.md). File layouts:
[file-formats.md](file-formats.md).

## Package and dependencies (P)

- Depends on Qt5 (core, gui, printsupport, concurrent), zlib, `r-cran-rcpp`; ships `RInside_0.2.13.tar.gz`
  (embedded R console), Qwt and QCustomPlot (plots), Boost (regex), ALGLIB (linear algebra).
- Ships one sample library `usr/share/metabolitedetector/lib_aa.lbr` (gzip, custom binary; amino acids).
- Maintainer `k.hiller@tu-bs.de`, TU Braunschweig.

## Workflow and features

| F | Feature | Old entry point (T/S) | Notes | New app |
|---|---------|----------------------|-------|---------|
| 1 | NetCDF (ANDI-MS) import | `&NetCDF import...`, `CDFBatchImporter`, `gcms::CdfImporter` | Batch converts `.cdf` to the internal `.bin`/`.idx` format; optional spectrum compression; optional TOF calibration. | done (read directly, no conversion; NetCDF export #30) |
| 2 | FastFlight import | `&FastFlight Import...`, `gcms::FastFlightImporter`, `MsData.FFT` | Raw TOF data; needs `.mzcalib` TOF calibration; can write `_centroid/_profile/_nominal.cdf`. | #29 |
| 3 | TOF m/z calibration | `TOFCalibrationWindow`, `TOFCalibrationCalculator` | Polynomial fit of calibration pairs from a reference spectrum; saved as `.mzcalib`. | #29 |
| 4 | Chromatogram viewer | `ChromatogramPlot`, `SpectrumPlot`, `MainWindow` | TIC, single-ion chromatogram, scan/time/ion spinners, slider, 1st-derivative plot, zoom, CSV/clipboard/image export, label compounds by score (≤0.7, >0.7, >0.9). | TIC/SIC, zoom, slider, unload done; derivative, time unit, info, score filter #79; export #9 |
| 5 | Spectrum view | `SpectrumPlot`, `SplitSpectrumPlot` | Profile/centroid, normalisation (none, base peak, unit, relative), split/difference mode for detected vs library spectrum. | normalisation #9; library comparison, centroid, component panel #77 |
| 6 | Chromatogram overlay | `ChromatogramOverlayDialog`, `MultiChromatogramPlot` | Overlay TIC or SIC of several runs. | done (#8) |
| 7 | Component detection (deconvolution) | `&Detect Components...`, `GCMSAnalyzer`, `PeakAnalyzer`, `DefaultPeakDetector`, `MatchedPeakDetector`, `ModelPeak`, `BaselineStatistics` | Peak detection per ion (derivative thresholds or matched filtering), baseline adjustment, noise units, binning (bins/scan), deconvolution width, min peaks, min % of base peak, model-peak integration. Results stored as `.cmp`. | default detector + deconvolution done; matched filtering #82; `.cmp` #13/#78 |
| 8 | Library identification | `gcms::LibrarySearch`, `SpectrumFunctions`, settings `LS_*` | Score = spectrum similarity, RI similarity or combined; pure/impure factors, cutoff, max ΔRI, mass filter ranges, scaled library, fragment pre-filter (`SQL_*`). | done (scaled library always on); `SQL_*` pre-filter and `LS_IDENTIFICATION_ALGORITHM` #21 |
| 9 | RI calibration wizard | `RICalibrationWizard`, `RICalculator` | Alkane (or other) calibration mix from a library, RT→RI table per chromatogram (`.calib`), optional "advanced" RT shift correction using a calibration compound; RT-only alignment as fallback. | basics + reference shift done (#20, #45); mix from library, editable table, per-run calibration, RT-only (#81) |
| 10 | Batch compound detection wizard | `Compound Detection Wizard...` | Runs #7 for many files after RI calibration. | targeted done (#22); per-run result list #78 |
| 11 | Batch quantification | `BatchQuantificationWizard`, `BatchQuantificator` | Targeted (reference library only) or non-targeted (temporary library collecting all spectra); replicate groups via regex; compound filter (reproducibility, discrepancy index, min ions); extended SIC rescan for missing quant ions; result table export. | targeted + groups done (#22); non-targeted #58; result views #83 |
| 12 | Quantification ions | `Compound::assignQuantificationIons`, `QuantIonChecker`, `QUANT_*` | N ions per compound, min distance, quality cut-off, excluded ions; manual/automatic ions in the library. | automatic ions in batch done; store in library #80, editor #12 |
| 13 | Normalisation | `CompoundSeriesNormalizer`, `PoolNormalizer` | By internal standard compound, summed sample signal, pool measurement, z-score. | #23 |
| 14 | Statistics | `BiomarkerDialog`, `tools::Statistics`, `BoxPlot`, `HistogramPlot` | Two-group biomarker table (means, ratio, p-value, hits), independent/paired samples, box plots. | #24 |
| 15 | PCA | `tools::PCA`, `PCAPlot` | Score plot with selectable PCs, loadings ("impact of metabolites"), cumulative variance. | #24 |
| 16 | MID wizard | `MIDWizard`, `tools::MIDSolver`, `NaturalMIDCalculator`, `MIDCompoundSeries` | Mass isotopomer distributions of labelled vs unlabelled runs; correction by sum formula, by unlabelled reference runs, by library spectrum, or uncorrected; M1 correction; tracer natural abundance; fragment syntax `(174/8)(261/5)`. | #28 |
| 17 | Natural isotope calculator | `NaturalIsotopeWidget`, `:/isotopes.csv` | Isotope pattern for a sum formula. | done (#27) |
| 18 | Library editor | `LibraryEditor`, `DialogAddLib`, `DialogImportLibrary`, `DialogExportLibrary` | Edit compounds, RI, quant ions, features (FRAG, comments); import/export AMDIS MSL (nominal/profile); add detected compound to library; copy RTs to a target library. | `.lbr` v2 read, MSL import done (#11), nominal MSL export in core only; export in the app and profile MSL #84; editor #12; add compound / copy RTs: core in `library_add.rs` (#80), app UI follows; v5 `.lbr` #41 |
| 19 | NIST search | `&NIST search...`, `gcms::NISTParser` | Search a user-supplied NIST library file. | #21 |
| 20 | R console | `RConsole`, `MDRConsole`, RInside | Results as R `data.frame data`, scripts, SVG/PDF plots, `.RData` export. | CSV + R template #26 (embedded R not ported on purpose) |
| 21 | Exports | `NetCDFExporter`, `CSVChromatogramExporter`, `CompoundExporter` | Chromatogram to NetCDF/CSV, compound tables to CSV, plots to PDF/PNG. | NetCDF done (#30); CSV/PNG/PDF/SVG/clipboard #9; compound CSV #78 |

## Parity with the new app

Audit of 2026-10-20 against `main` (#76). "done" means the new app covers it; `#N` is the open issue
that owns it. Sub-features below come from the binary's menus, dialogs and settings; each one is owned
by exactly one issue.

| F | Sub-feature | New app |
|---|-------------|---------|
| 1 | Batch conversion, spectrum compression, recursive search | not needed: runs are read directly |
| 2, 3 | FastFlight `.FFT`/`.FFC` import, centroiding (`SF_*`), TOF calibration window, `.mzcalib`, Calibration Setter | #29 |
| 4 | First derivative, time unit (ms/s/min), run info dialog, label filter by score, auto-scale Y | #79 |
| 4, 21 | Copy data/image, CSV, PNG/PDF/SVG, image size dialog | #9 |
| 5 | Split/difference/normal comparison with the library spectrum, centroid toggle, spectrum export | #77 |
| 5 | Component panel: S/N, sharpness, discrepancy index, scores, ion counts, peak list | #77 |
| 7 | Matched-filter peak detection for components | #82 |
| 7 | Per-run result list: table, find compound, CSV, save/load `.cmp`, `.cmp` as batch input | #78 |
| 8, 19 | NIST import into SQLite, fragment pre-filter (`SQL_*`), `LS_IDENTIFICATION_ALGORITHM`, NIST path | #21 |
| 11 | Non-targeted batch, rescan for missing quant ions, excluded ions, discrepancy index, S/N and required-peaks filters | #58 |
| 11 | Result tree, alternatives, RT table, rel. std. error, blanks / P(blank), delete rows/columns, all-ions CSV | #83 |
| 12, 18 | Store detected compound and quant ions in the library, copy RTs to a target library | #80 |
| 18 | Editor: auto quant ions, features dialog, ion frequency, resolution, MID fragments | #12 |
| 18 | MSL export in the app, also as profile/FWHM spectrum | #84 |
| 18 | `.lbr` version 5 | #41 |
| 13 | Internal standard, summed signal, pool (with pool mapping export), z-score | #23 |
| 14, 15 | Group assignment, p-value filter, box plots, PCA axes/loadings/variance, plot export | #24 |
| 16 | MID wizard: methods, tracer, M1, result tree, StdErr, CSV/clipboard, MID plots | #28 |
| 20 | CSV results + R template instead of the embedded R console | #26 |
| files | `.cmp` and binary `.calib` layout, RT unit | #13 |
| files | High-resolution runs (`getDataType() == 2`) | #34 |

Issues that waited for lab or licensed files (#13, #21, #28, #29, #34, #41) now have a plan on the issue to produce references with the
old binary through the golden harness (`tools/golden/`) or with synthetic data, so none of them waits
for lab files any more.

## Internal file types (T)

`.bin` + `.idx` converted chromatogram · `.cmp` detected compounds · `.lbr` compound library (gzip) ·
`.calib` RI calibration · `.mzcalib` TOF calibration · `.msl` AMDIS library (import/export) ·
NIST library (import for search) · settings stored by `GCMSSettings::toDisk`.

## Observations about the old code quality (T/S)

- Long-running work on Qt threads with ad-hoc locks (`CRITICAL: Was not able to obtain lock in MainWindow::setLibrary()`).
- Mixed German/English UI text (`Keine Lib`), hard-coded home paths (`/home/kh`).
- R integration via RInside makes Linux packaging depend on a matching R installation – the main
  reason the app does not run on Windows.
