# Batch quantification: untargeted mode, filters and the rescan of missing ions (issues #58, #110)

Source: the old Linux binary 3.3 (`/usr/bin/metabolitedetector_app`), read with objdump in WSL because the
`rea-wsl` decompiler crashed on these functions. Addresses are `nm` addresses of that binary (Ghidra = nm +
0x100000). "Observed" is read from the disassembly; "Not determinable" is what the disassembly does not settle.
The port is `crates/md-core/src/batch.rs` and its submodules `batch/untargeted.rs` and `batch/rescan.rs`;
`md-cli` and the app only call them (`complete_scan`). The golden harness (`tools/golden/`) drives the single-run
pipeline only; there is no way to feed it a batch, so nothing here is compared against live output of the old
program (Observed from the code only).

## Pipeline in the old program (`BatchQuantificator::startAnalysis`, 0x472058)

1. Per run, `createResultSeries` / `createReplicateResultSeries` (0x6358be, 0x63b4c2) load the run's component
   file and apply the component filter (below) to its components, then match them against the library
   (`addSeriesFirstStage`). Targeted: the library is the reference library. Untargeted
   (`BatchQuantificator+0x28` is null): a temporary library made by `createLibrary` (0x47334e).
2. `setReplicateGroupMapping`, then `removeNonReproducableCompounds(reproducibility)` (0x63af4c).
3. Only if "Extended SIC scan" is ticked (QSettings `bsw_completion`): for every run `GCMSDiskScan` +
   `PeakAnalyzer` + `ResultSeries::preciseSeriesCompletion` (0x635006), the rescan.

So the order is filter, match, reproducibility, rescan. The rescan comes last and nothing filters its values.

## Settings of the dialog (`BatchQuantificatorSettings`, stored at `BatchQuantificator+0x38`; Observed)

| Field | Widget / key | Default |
|---|---|---|
| ΔRI (`LS_RI_DIFF`) | | 20 |
| score cutoff of the library search | "Req. Score" | |
| Compound reproducibility | `bsw_rep` | 0.5 |
| Req. S/N | `bsw_sn` | 0.0 |
| Max. peak discr. index | `bsw_dis` | 100 |
| Min # ions | `bsw_req_peaks` | 1 |
| Extended SIC scan | `bsw_completion` | |
| Redetect all quant. ions | `bsw_redetect` | |

"Mass filter in use" is the `QUANT_EXCL_IONS` list being non-empty.

## Component filter (`ResultSeries::filterCompounds` 0x636fb0, `BatchQuantificator::filterCompounds` 0x475876; Observed)

- Arguments `(components, discrepancy limit, min ions, S/N limit)`. A component is kept iff
  `ions > min_ions`, `mean peak quality < limit` (strict) and `mean S/N > limit` (strict) or the mean S/N is
  NaN. The means are over all ions of the component (`getAveragePeakQuality`, `getAverageSignalToNoise`).
- It applies **per component of each run**, not per compound: the load of the run drops the failing ones
  before they are matched. A run whose component fails is therefore a run without that compound, and the
  rescan may fill it afterwards.
- With the dialog's filter switched off the values are `1`, `1000.0` and `0.0`.
- The copy used while building the temporary library keeps `ions >= min_ions` and `S/N >= limit` (the
  boundaries differ by one step from the series version) and prints "Erased a compound".

## `preciseSeriesCompletion` (Observed)

- The series is looked up by name; an unknown name throws `CompoundSeriesException("No such name")`.
- Skipped (`INFORMATION: Skipping this series, because only N samples have an integral!=0.0`) when fewer than
  `(long)(n / 4.0 + 0.5)` samples have a summed integral != 0, `n` being the sample count (constants 4.0 at
  0xa96d68, 0.5 at 0xa96c50). That is about a quarter of the runs, not half.
- For every sample with a zero integral and every quant ion whose `getIntegral(ion, sample)` is 0 it calls
  `scanIons`. If the result differs from 0 it prints `INFORMATION: Found additional information for compound X
  Ion: Y Index: i` and stores it with `CompoundSeries::setIntegral`.
- **Window.** The first argument is ΔRI of the dialog (`startAnalysis` passes `BatchQuantificator+0x38`).
  - RI branch: used when the series has an average RI >= 0 and `RICalculator::isCalibrated()` of the scanned
    run. Positions are the series' RIs >= 0, width = ΔRI, last argument 0.
  - RT branch: positions are the series' retention times >= 0 (ms), width = ΔRI * 1000.0 (0xa96ce8), last
    argument 1. So the window is ΔRI seconds wide in total.
  - The window is `min(positions) - width/2 .. max(positions) + width/2`.

## `ResultSeries::scanIons` (0x637e30..0x639018, Observed)

- Takes the chromatogram of the ion (`GCMSDiskScan::getChromatogram`), runs a fresh `DefaultPeakDetector`, keeps
  the peaks whose ion lies within +-0.5 of the wanted one and whose apex scan lies in the window.
  `INFORMATION: No peak found for ion: X` if none.
- Builds a `GCMSAnalyzer` on the kept peaks with `setMinimumPeakNumberRequired(1)` and calls
  `deconvoluteChromatogram(peaks, scan_count, start_scan, 2 * AN_DECONVOLUTION_WIDTH)`, once per ion.
- **Peak choice (#336).** The expected position is the **mean** of the series' positions (the same vector the
  window comes from), not the middle of the window. `std::sort` with `Compound_Compare(mean, ion)` orders the
  components by `|RI - mean|` when both have an RI >= 0, else by `|RT - mean|`; on equal distance the one with
  the larger intensity of the ion (`getIntensity(ion, 0.5)`) comes first (`Compound_Compare::operator()`,
  0x6b9c34). The 0.5 is the m/z tolerance of that lookup, **not** a limit on the distance: any component in
  the window can be taken, however far from the mean.
- It walks the sorted components and returns the first non-zero integral of the ion (`getIntensity(ion, 0.5)`
  of the component's integral spectrum); a component without the ion prints a message and the next is tried.
- If no component has the ion, the peaks are grouped per ion value; per group it keeps the one peak whose apex
  scan is nearest to `rt2Scan(mean)` (`ri2Scan` in the RI branch, first one on a tie), sums the raw values of
  the kept peaks (`Peak::getValues`) and prints `INFORMATION: Found additional peak areas without deconvolution`.
- Port (`batch/rescan.rs`): per missing ion the peak nearest to the mean, the higher one on a tie. So with
  ΔRI 20 and no RI calibration, Glycine 3TMS (313 s) gets the m/z 177 peak at 321 s in the example runs, as in
  3.3: its own `DefaultPeakDetector` finds no other m/z 177 peak in 302.9..323.1 s (`old-md-3.3.tsv`, runs
  7235eg04, 20 and 32). Left as is: the port deconvolutes the picked peaks of all missing ions together and
  does not fall back to the next component when the fit loses the ion (it sums the raw values of the peak).
- `ResultSeriesException("Not implemented yet")` is thrown when the byte at `ResultSeries+0x141` is false. The
  constructor stores its first bool there and both `createResultSeries` and `createReplicateResultSeries` pass 1,
  so batch quantification never reaches it.

## `removeNonReproducableCompounds` (0x63af4c, Observed)

A series stays if some replicate group has `getGroupHits / getGroupMemberCount >= reproducibility` (the loop
leaves at the first group that passes); otherwise `removeCompoundSeries`. Runs before the rescan.

## Untargeted: temporary library (`createLibrary`, `addCompoundSpectra`, Observed)

- For each replicate group a new `LibrarySearch`; for each file of the group its components are filtered and
  passed to `addCompoundSpectra(components, cutoff, source = file name)`.
- `addCompoundSpectra` (0x6280f0) calls `mapCompoundsToLibrary` (0x629526; one-to-one, best score first, see
  algorithms.md §5.7) for **all** components of the file against the library as it is at that moment. For each
  match with an overall score >= cutoff and a source not yet in the library entry, the spectrum is merged into it
  (`LibraryCompound::addCompoundSpectrum`, 0x61cf4a) and removed from the set; a spectrum from a source already
  contained is skipped ("Not adding this compound spectrum since it was already contained from this source").
  What is left becomes new library entries, named `name_N` when the name exists.
- Merge: per ion the mean over the spectra that contain the ion, `calculateNewValue(nA, nB, vA, vB) =
  nA/(nA+nB) * vA + nB/(nA+nB) * vB`, cut to an integer for the `<int,float>` library (the +0.5 applies only to
  float-valued spectra). New ions are appended. Retention time and RI are weighted by the spectrum count.
- After all files of a group, entries with `spectrumCount / filesInGroup < reproducibility` are removed. The
  group libraries are then merged into the global one with the same call (set of entries instead of
  components), the quant ions are set from `getCentroidLibraryCompound`, and entries are named from an optional
  reference library (features `lib_score`, `source_file`, `CASNO`, `best_lib_score`, `best_matching`,
  `alt_ref_hits`; names like "No direct match: ...").
- **Not determinable:** how the mass tolerance 0.001 (0xa96d58) and an RI of -1 enter the merge. The port
  works on nominal masses, where the tolerance changes nothing; for the RI it lets the known value win.

## Retention time shift

The old batch used the RT shift correction (§6.5) per run like the single-run search. The port does the same:
`ScanSettings::shift` shifts the calibration onto the reference substance in each run, `RunResult::calibration`
keeps the shifted calibration, and the rescan converts the expected positions between runs through the RI.

## ANOVA p-value of the group table (`Statistics::oneWayANOVA` 0x6c3af6, `getANOVAPValue` 0x63ac2c; Observed)

The "p-value" column of the old group statistics is a one-way ANOVA over the replicate groups of a substance, blank
group included. `getANOVAPValue` returns 1 for one group or none; otherwise `oneWayANOVA` takes the values of
each group (missing = 0, normalised if a normalisation is set) and computes with `N` values in `k` groups:

```
m = mean(all);   TSS = Σ (x - m)²;   SSB = Σ nᵢ (meanᵢ - m)²;   SSW = TSS - SSB
F = (SSB / (k - 1)) / (SSW / (N - k));   p = gsl_cdf_fdist_Q(F, k - 1, N - k)   // upper tail
```

The port: `stats::anova` (tail via the regularized incomplete beta function), shown as the last column of
"Group statistics" (`P (ANOVA)`). `None` ("–") where the old code would divide by zero: an empty group, `N = k`,
no spread within or between groups.

Verified (#222): `tools/anova/old_anova.cpp` (`make-anova.sh`, `LD_PRELOAD` like the MID tools) calls `oneWayANOVA`
on 74 cases (two of them never reach it in `getANOVAPValue`; hand cases, unequal groups, zeros, p down to 6e-23, 60 random); `crates/md-core/tests/old_anova.rs`
matches the other 72 to 1e-9 relative (checked again with the final generator). The old function divides 0 by 0 for identical values (nan) and, with
`N = k` or an empty group, ends in the GSL error handler (the whole program aborts); the port gives no p-value
("–") in these cases. `getANOVAPValue` never calls it with fewer than two groups.

Export: `BatchQuantificationResultWidget::exportData` (0x4506c0) writes the shown table model cell by cell (header
row, then the rows of the model), so the group view with its p-value column is exported as shown. The port writes
`FILE.groupstats.csv` next to the results CSV (`Table::write_group_stats_file`).

## Port: what differs

| Old | Port |
|---|---|
| One temporary library per replicate group, merged afterwards, reproducibility inside | One pooled library over all runs (`batch/untargeted.rs`); groups and reproducibility are applied by `table`. Within a run the snapshot rule, one component per entry and the mean spectrum are the same. Entries are named `Unknown n`. |
| RI of -1 in the merge | The known RI wins (`shortcut:` in `untargeted.rs`). |
| Component filter at load, before matching and the rescan | Applied by `table` on the ions of the quant ions (S/N per ion, mean quality of the quant ions, `>=`), so it can be changed without a new scan. Differences: a run filtered out is not rescanned, and rescanned values are filtered too. Defaults: off. |
| Reproducibility before the rescan | `table` counts a run only if its value exists without the ions the rescan added (`Found::rescued`, `Found::rescanned`). Whether `getGroupHits` counts identifications or integrals is not determinable from the disassembly; the port counts runs with a value. |
| Rescan window ΔRI wide in RI (calibrated run) or ΔRI s wide in RT | Same (ΔRI is `SearchSettings::ri_diff`); the positions of the other runs are the ones with an RI under their own calibration, else their retention times. |
| Rescan of the ions in the order of the series | Quant ions as `table` chooses them with the same settings; changing them afterwards needs a new rescan. |
| Rescanned substance absent in a run | `Found { rescanned: true, score: 0 }`: counts for values and reproducibility, not for the mean retention time or the choice of quant ions. |
