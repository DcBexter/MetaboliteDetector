# Roadmap

Milestones for the Rust rewrite. Each bullet is a GitHub issue; the [project board](https://github.com/users/DcBexter/projects/2)
owns status and order, the issues own scope, acceptance and dependencies. Feature numbers (F1–F21) refer to
[reverse-engineering/feature-inventory.md](reverse-engineering/feature-inventory.md).

Every milestone keeps the [performance rules](architecture.md#performance-rules): worker threads,
parallel loops, pixel-bounded drawing and a perf case per new algorithm.

## M1 – Viewer

[Milestone](https://github.com/DcBexter/MetaboliteDetector/milestone/2)

- [Workspace, NetCDF import, TIC/spectrum viewer, example data (F1, F4)](https://github.com/DcBexter/MetaboliteDetector/issues/7)
- [SIC by m/z, overlay of several runs (F4, F6)](https://github.com/DcBexter/MetaboliteDetector/issues/8)
- [Spectrum normalisation, CSV/PNG export (F5, F21)](https://github.com/DcBexter/MetaboliteDetector/issues/9)
- [CI on Windows and Linux](https://github.com/DcBexter/MetaboliteDetector/issues/10)

## M2 – Libraries and compatibility

[Milestone](https://github.com/DcBexter/MetaboliteDetector/milestone/3)

- [Read `.lbr`, import/export MSL (F18)](https://github.com/DcBexter/MetaboliteDetector/issues/11)
- [Library editor (F18)](https://github.com/DcBexter/MetaboliteDetector/issues/12)
- [Read old `.cmp` and `.calib`](https://github.com/DcBexter/MetaboliteDetector/issues/13)

## M3 – Component detection (core science)

[Milestone](https://github.com/DcBexter/MetaboliteDetector/milestone/4)

- [Reverse-engineer deconvolution, scoring, RI, quant ions](https://github.com/DcBexter/MetaboliteDetector/issues/14)
- [Baseline, noise and peak detection (F7)](https://github.com/DcBexter/MetaboliteDetector/issues/15)
- [Deconvolution and integration (F7)](https://github.com/DcBexter/MetaboliteDetector/issues/16)
- [Golden tests against the old app](https://github.com/DcBexter/MetaboliteDetector/issues/17)
- [Settings with the old parameter set](https://github.com/DcBexter/MetaboliteDetector/issues/18)

## M4 – Identification and RI

[Milestone](https://github.com/DcBexter/MetaboliteDetector/milestone/5)

- [Library search (F8)](https://github.com/DcBexter/MetaboliteDetector/issues/19)
- [RI calibration wizard (F9)](https://github.com/DcBexter/MetaboliteDetector/issues/20)
- [NIST search (F19)](https://github.com/DcBexter/MetaboliteDetector/issues/21)

## M5 – Batch quantification and statistics

[Milestone](https://github.com/DcBexter/MetaboliteDetector/milestone/6)

- [Batch detection and quantification (F10–F12)](https://github.com/DcBexter/MetaboliteDetector/issues/22)
- [Normalisation (F13)](https://github.com/DcBexter/MetaboliteDetector/issues/23)
- [Biomarker statistics, box plots, PCA (F14, F15)](https://github.com/DcBexter/MetaboliteDetector/issues/24)
- [Headless CLI](https://github.com/DcBexter/MetaboliteDetector/issues/25)
- [CSV/R export replacing the R console (F20)](https://github.com/DcBexter/MetaboliteDetector/issues/26)

## M6 – Isotope labelling

[Milestone](https://github.com/DcBexter/MetaboliteDetector/milestone/7)

- [Natural isotope calculator (F17)](https://github.com/DcBexter/MetaboliteDetector/issues/27)
- [MID wizard (F16)](https://github.com/DcBexter/MetaboliteDetector/issues/28)

## Later / low priority

[Milestone](https://github.com/DcBexter/MetaboliteDetector/milestone/8)

- [FastFlight TOF import and m/z calibration (F2, F3)](https://github.com/DcBexter/MetaboliteDetector/issues/29)
- [Old `.bin`/`.idx` reader, NetCDF export](https://github.com/DcBexter/MetaboliteDetector/issues/30)

## Reverse-engineering status

| Area | Status | Document |
|------|--------|----------|
| Feature inventory | done | [feature-inventory.md](reverse-engineering/feature-inventory.md) |
| Noise, baseline, peak detection | done (pseudocode, not yet validated against the old app) | [algorithms.md](reverse-engineering/algorithms.md) §2–3 |
| Deconvolution, library scoring, RI calibration, quantification ion choice | done; open details listed per section | [algorithms.md](reverse-engineering/algorithms.md) §4–7 |
| MID, normalisation, statistics, PCA | open | [algorithms.md](reverse-engineering/algorithms.md) §8 |
| File formats and settings | mostly done: `.lbr` v2 parser verified on the sample library; `.bin`, `.cmp`, `.calib`, `.mzcalib`, settings layouts documented. Open: `.idx` layout, peak records in `.cmp`, MSL export | [file-formats.md](reverse-engineering/file-formats.md) |
