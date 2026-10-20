# Batch quantification: untargeted mode, filters and the rescan of missing ions (issue #58)

Source: the old Linux binary 3.3 (`/usr/bin/metabolitedetector_app`), read with objdump in WSL because the
`rea-wsl` decompiler crashed on these functions. Addresses are virtual addresses of that binary. "Observed" is
read from the disassembly; "Not decoded" is a gap that the port fills with an own choice, marked
`shortcut:` in the code. The port is `crates/md-core/src/batch.rs` and its submodules `batch/untargeted.rs`
and `batch/rescan.rs`; `md-cli` and the app only call them (`complete_scan`).

## Pipeline in the old program (Observed)

`BatchQuantificator` scans each run (peaks, deconvolution, `LibrarySearch`), collects the results per compound
in a `ResultSeries` (one `CompoundSeries` per compound, one entry per sample) and then:

1. Targeted: the compounds are the entries of the reference library. Non-targeted: a temporary library
   (`BatchQuantificator::createLibrary`, `LibrarySearch::addCompoundSpectra`) collects the spectra of all
   components of all runs, and every run is then matched against it.
2. `ResultSeries::preciseSeriesCompletion` (0x635006) fills the missing ions (below).
3. `BatchQuantificator::filterCompounds` drops compounds that fail the filters.

## `preciseSeriesCompletion` (Observed)

- The series is looked up by name; an unknown name throws `CompoundSeriesException("No such name")`.
- A series is skipped with `INFORMATION: Skipping this series, because only N samples have an integral!=0.0`
  when too few samples have a summed integral. The threshold is computed from the sample count and a constant
  0.5 next to a 4.0 in rodata (0xa96c50, 0xa96d68). The port uses "found in at least half of the runs"
  (`RescanSettings::min_share`); the exact comparison is **Not decoded**.
- For every sample with a zero integral (`CompoundSeries::hasZeroIntegrals`) and every quant ion
  (`getQuantificationIons`) whose `getIntegral(ion, sample)` is 0, the search window is built from the series:
  the retention indices when the average RI is known and the run is calibrated
  (`RICalculator::isCalibrated`), else the retention times. Only values >= 0 are used.
- `ResultSeries::scanIons(analyzer, positions, ion, width, use_ri)` returns the integral. If it differs from
  0, the program prints `INFORMATION: Found additional information for compound X Ion: Y Index: i` and stores it
  with `CompoundSeries::setIntegral(ion, sample, value)`.

## `ResultSeries::scanIons` (0x637e30..0x639018, Observed unless marked)

- Window: `min(positions) - width/2` .. `max(positions) + width/2` (division by the constant 2.0); the flag
  chooses `rt2Scan` or `ri2Scan` to turn it into scans. The mean of the positions is the expected position.
- Takes the chromatogram of the ion (`GCMSDiskScan::getChromatogram`), runs a fresh `DefaultPeakDetector`, keeps
  the peaks whose ion lies within +-0.5 of the wanted one and whose apex scan lies in the window.
  `INFORMATION: No peak found for ion: X` if none.
- Builds a `GCMSAnalyzer` on the kept peaks with `setMinimumPeakNumberRequired(1)` and calls
  `deconvoluteChromatogram(peaks, scan_count, start_scan, 2 * AN_DECONVOLUTION_WIDTH)` (0x638589..0x6385d7).
  From the components it takes the one closest to the expected position (`Compound_Compare(double, float)`,
  tolerance 0.5) and returns the integral of the ion in it.
- If the ion is not in that component, the peaks are grouped per ion (`map<float, list<Peak*>>`), the raw
  values of the peaks in the window are summed (`Peak::getValues`) and
  `INFORMATION: Found additional peak areas without deconvolution` is printed. One more branch throws
  `ResultSeriesException("Not implemented yet")`; which input reaches it is **Not decoded**.
- **Not decoded:** the source of the `width` argument (a constant 1000.0 and a stored double appear next to the
  calls) and so its unit. The port searches the expected times of the other runs +- `RescanSettings::window_s`/2
  (default 6 s in total) in the run's retention time.

## Filters (`BatchQuantificator::filterCompounds`)

- Reproducibility: a compound stays if some replicate group has `hits / members >= threshold`.
- Minimum number of quant ions found per sample.
- Discrepancy index: mean peak quality of the ions of the compound (the "Discrepancy" of the component panel,
  §3.7 `determinePeakQuality`).
- Required S/N and required peaks (the latter is `AN_MIN_PEAK_NUMBER` of the scan).

**Not decoded:** the argument order of `ResultSeries::filterCompounds`, whether the discrepancy and S/N limits apply per
sample or per compound, the default of the S/N limit, and whether the filter runs before or after
`preciseSeriesCompletion`. The port applies S/N and discrepancy per cell (an ion counts only above the S/N limit;
a cell counts only if the mean quality of its ions is at most the limit), after the rescan, as part of `table`.
Both limits are off by default.

## Retention time shift

The old batch used the RT shift correction (§6.5) per run like the single-run search. The port does the same:
`ScanSettings::shift` shifts the calibration onto the reference substance in each run, `RunResult::calibration`
keeps the shifted calibration, and the rescan converts the expected times between runs through the RI.

## Port: what differs

| Old | Port |
|---|---|
| Temporary library built by `addCompoundSpectra` | Greedy pooling in run order: a component joins the group with the best library search score above the cutoff (spectrum, plus RI when calibrated), at most one component per run and group; a group keeps its first spectrum (`batch/untargeted.rs`). Groups are named `Unknown n`. |
| Rescan per ion with window from the series | Same idea; one local deconvolution per substance and run with the peaks of its missing ions (min. 1 peak, 2 x width); the raw peak sum when the fit gives nothing. |
| Rescan of the ions in the order of the series | Quant ions as `table` chooses them with the same settings; changing them afterwards needs a new rescan. |
| Rescanned substance absent in a run | `Found { rescanned: true, score: 0 }`: counts for values and reproducibility, not for the mean retention time or the choice of quant ions. |
