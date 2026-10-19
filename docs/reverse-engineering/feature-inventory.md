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

| F | Feature | Old entry point (T/S) | Notes |
|---|---------|----------------------|-------|
| 1 | NetCDF (ANDI-MS) import | `&NetCDF import...`, `CDFBatchImporter`, `gcms::CdfImporter` | Batch converts `.cdf` to the internal `.bin`/`.idx` format; optional spectrum compression; optional TOF calibration. |
| 2 | FastFlight import | `&FastFlight Import...`, `gcms::FastFlightImporter`, `MsData.FFT` | Raw TOF data; needs `.mzcalib` TOF calibration; can write `_centroid/_profile/_nominal.cdf`. |
| 3 | TOF m/z calibration | `TOFCalibrationWindow`, `TOFCalibrationCalculator` | Polynomial fit of calibration pairs from a reference spectrum; saved as `.mzcalib`. |
| 4 | Chromatogram viewer | `ChromatogramPlot`, `SpectrumPlot`, `MainWindow` | TIC, single-ion chromatogram, scan/time/ion spinners, slider, 1st-derivative plot, zoom, CSV/clipboard/image export, label compounds by score (≤0.7, >0.7, >0.9). |
| 5 | Spectrum view | `SpectrumPlot`, `SplitSpectrumPlot` | Profile/centroid, normalisation (none, base peak, unit, relative), split/difference mode for detected vs library spectrum. |
| 6 | Chromatogram overlay | `ChromatogramOverlayDialog`, `MultiChromatogramPlot` | Overlay TIC or SIC of several runs. |
| 7 | Component detection (deconvolution) | `&Detect Components...`, `GCMSAnalyzer`, `PeakAnalyzer`, `DefaultPeakDetector`, `MatchedPeakDetector`, `ModelPeak`, `BaselineStatistics` | Peak detection per ion (derivative thresholds or matched filtering), baseline adjustment, noise units, binning (bins/scan), deconvolution width, min peaks, min % of base peak, model-peak integration. Results stored as `.cmp`. |
| 8 | Library identification | `gcms::LibrarySearch`, `SpectrumFunctions`, settings `LS_*` | Score = spectrum similarity, RI similarity or combined; pure/impure factors, cutoff, max ΔRI, mass filter ranges, scaled library, fragment pre-filter (`SQL_*`). |
| 9 | RI calibration wizard | `RICalibrationWizard`, `RICalculator` | Alkane (or other) calibration mix from a library, RT→RI table per chromatogram (`.calib`), optional "advanced" RT shift correction using a calibration compound; RT-only alignment as fallback. |
| 10 | Batch compound detection wizard | `Compound Detection Wizard...` | Runs #7 for many files after RI calibration. |
| 11 | Batch quantification | `BatchQuantificationWizard`, `BatchQuantificator` | Targeted (reference library only) or non-targeted (temporary library collecting all spectra); replicate groups via regex; compound filter (reproducibility, discrepancy index, min ions); extended SIC rescan for missing quant ions; result table export. |
| 12 | Quantification ions | `Compound::assignQuantificationIons`, `QuantIonChecker`, `QUANT_*` | N ions per compound, min distance, quality cut-off, excluded ions; manual/automatic ions in the library. |
| 13 | Normalisation | `CompoundSeriesNormalizer`, `PoolNormalizer` | By internal standard compound, summed sample signal, pool measurement, z-score. |
| 14 | Statistics | `BiomarkerDialog`, `tools::Statistics`, `BoxPlot`, `HistogramPlot` | Two-group biomarker table (means, ratio, p-value, hits), independent/paired samples, box plots. |
| 15 | PCA | `tools::PCA`, `PCAPlot` | Score plot with selectable PCs, loadings ("impact of metabolites"), cumulative variance. |
| 16 | MID wizard | `MIDWizard`, `tools::MIDSolver`, `NaturalMIDCalculator`, `MIDCompoundSeries` | Mass isotopomer distributions of labelled vs unlabelled runs; correction by sum formula, by unlabelled reference runs, by library spectrum, or uncorrected; M1 correction; tracer natural abundance; fragment syntax `(174/8)(261/5)`. |
| 17 | Natural isotope calculator | `NaturalIsotopeWidget`, `:/isotopes.csv` | Isotope pattern for a sum formula. |
| 18 | Library editor | `LibraryEditor`, `DialogAddLib`, `DialogImportLibrary`, `DialogExportLibrary` | Edit compounds, RI, quant ions, features (FRAG, comments); import/export AMDIS MSL (nominal/profile); add detected compound to library; copy RTs to a target library. |
| 19 | NIST search | `&NIST search...`, `gcms::NISTParser` | Search a user-supplied NIST library file. |
| 20 | R console | `RConsole`, `MDRConsole`, RInside | Results as R `data.frame data`, scripts, SVG/PDF plots, `.RData` export. |
| 21 | Exports | `NetCDFExporter`, `CSVChromatogramExporter`, `CompoundExporter` | Chromatogram to NetCDF/CSV, compound tables to CSV, plots to PDF/PNG. |

## Internal file types (T)

`.bin` + `.idx` converted chromatogram · `.cmp` detected compounds · `.lbr` compound library (gzip) ·
`.calib` RI calibration · `.mzcalib` TOF calibration · `.msl` AMDIS library (import/export) ·
NIST library (import for search) · settings stored by `GCMSSettings::toDisk`.

## Observations about the old code quality (T/S)

- Long-running work on Qt threads with ad-hoc locks (`CRITICAL: Was not able to obtain lock in MainWindow::setLibrary()`).
- Mixed German/English UI text (`Keine Lib`), hard-coded home paths (`/home/kh`).
- R integration via RInside makes Linux packaging depend on a matching R installation – the main
  reason the app does not run on Windows.
