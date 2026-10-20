# MetaboliteDetector 3.3 (Linux x86-64) - algorithm notes

Status: partial. Noise/baseline, peak detection and peak metrics (ported in `crates/md-core/src/peaks.rs`, deviations in `docs/architecture.md#peak-detection`), deconvolution (Section 4), library scoring (Section 5, ported in `crates/md-core/src/search.rs`), RI calibration (Section 6) and quantification ion selection (Section 7) are reconstructed. Normalisation (Section 7c) is reconstructed. Statistics and PCA (Section 7f) are reconstructed. MID (Section 7g) is reconstructed.

Conventions
- Observed = read from decompiler or disassembly output. Inferred = deduced from usage, names or typical patterns; verify before porting.
- Addresses: objdump / symbol-table (`syms`) addresses are file addresses. Ghidra MCP addresses = syms address + 0x100000. Example: `calculateBaseline` is syms 0x5ed5ee, MCP 0x6ed5ee.
- Template instances: `<int,float>` (int intensities, float derivative/weights). Other instances (`<double,double>`) exist; not compared.
- Ghidra caveat: decompiles of functions that call `operator new` are truncated after that call (no-return thunk). Control flow for those was recovered from objdump.
- Ghidra local variable names are unreliable (off-by-one stack slots). Names below are from usage.

## 1. Settings persistence (Observed) - needed to understand "defaults"

- `gcms::GCMSSettings::toStream(std::ostream&)` (syms 0x6d5c1c), `fromStream(std::istream&)` (0x6d6052), `toDisk(char const*)` (0x6d64c6).
- The file is written and read through a gzip stream wrapper (`gzstreambase` / `gzstreambuf`, gzopen/gzread/gzwrite).
- Format: text. For each key, toStream writes the key name, then the value; fromStream reads with `>>` in the same order. Bools are written as single bytes via movzbl/mov %al.
- Key order (as written by toStream, Observed):
  1. SF_PEAK_THRESHOLD_BEGIN, SF_PEAK_THRESHOLD_END, SF_MAX_BASE_DEV, SF_FWHM
  2. SF_CENTROID_THRESHOLD_BEGIN, SF_CENTROID_THRESHOLD_END, SF_CENTROID_MAX_BASE_DEV
  3. AN_PEAK_THRESHOLD_BEGIN, AN_PEAK_THRESHOLD_END, AN_MAX_BASE_DEV, AN_MIN_PEAK_HEIGHT
  4. AN_DECONVOLUTION_WIDTH, AN_NUMBER_BINS, AN_MIN_PEAK_NUMBER
  5. LS_RI_DIFF, LS_PURE_FACTOR, LS_IMPURE_FACTOR, LS_SCALE_FACTOR_ION, LS_SCALE_FACTOR_VALUE, LS_SCALE_FACTOR_NO_MATCH, LS_RESOLUTION
  6. LS_USE_RI (bool), LS_USE_SCALED_LIB (bool), LS_USE_SPEC (bool)
  7. QUANT_ION_NUMBER, QUANT_MIN_DIST_IONS, QUANT_QUALITY_CUTOFF
  8. QUANT_EXCL_IONS (std::set<pair<float,float>>; written as size then pairs)
- Keys also seen in code but not in toStream: AN_MIN_HEIGHT_BP, AN_MODELED_VALUES_FOR_INTEGRATION, SQL_MINIMAL_NUMBER_EQUAL_FRAGMENTS, SQL_MAX_NUMBER_FRAGMENTS, LS_CUTOFF, LS_MASS_FILTER, LS_IDENTIFICATION_ALGORITHM. The SQL_* keys are QSettings values and
LS_IDENTIFICATION_ALGORITHM is never written (file-formats.md §5); the others' persistence is Unknown.
- Default values (Observed): the static members are initialised in `.data`. The values are AN_PEAK_THRESHOLD_BEGIN 10, AN_PEAK_THRESHOLD_END -5, AN_MAX_BASE_DEV 20, AN_MIN_PEAK_HEIGHT 20, AN_NUMBER_BINS 10, AN_DECONVOLUTION_WIDTH 1.0, AN_MIN_PEAK_NUMBER 10 and AN_MODELED_VALUES_FOR_INTEGRATION true. AN_MIN_HEIGHT_BP is in `.bss` (0). A fresh `~/.config/TU-Braunschweig/MetaboliteDetector.conf` holds the same values. The GCMSSettings ctor (syms 0x3dfb98) is empty.

Settings used by the analysis (name, storage address, where read):
- AN_PEAK_THRESHOLD_BEGIN (0xcc20e0), AN_PEAK_THRESHOLD_END (0xcc20e8): slope thresholds in analyzeBaseline and DefaultPeakDetector.
- AN_MIN_PEAK_HEIGHT (0xcc20f8): S/N threshold in calculateBaseline.
- AN_NUMBER_BINS (0xcc2100): time quantization in determinePreciseMaximizationTime.
- AN_DECONVOLUTION_WIDTH (0xcc2108), AN_MIN_PEAK_NUMBER (0xcc2110): used by deconvolution (not analyzed).
- SF_PEAK_THRESHOLD_BEGIN (0xcc20a8), SF_PEAK_THRESHOLD_END (0xcc20b0), SF_FWHM (0xcc20c0), SF_MAX_BASE_DEV (0xcc20b8): spectrum-filter settings (usage not analyzed).
- AN_MAX_BASE_DEV (0xcc20f0): purpose Unknown.

## 2. Noise and baseline (Observed unless marked)

### 2.1 analyzeBaseline<int,float> (MCP 0x6ee06e; syms 0x5ee06e)

Called via a virtual slot (vtbl[5]) from DefaultPeakDetector::detectPeaks with (SG-smoothed intensities, derivative, segment_length = 15).

```
analyzeBaseline(values v[N], slope s[N], L):
  segments: start at 0, step L, while i < N - L (non-overlapping)
  accept segment if for every k in segment:
      AN_PEAK_THRESHOLD_END < s[k] < AN_PEAK_THRESHOLD_BEGIN
  for accepted segment:
      m   = mean(v over segment)
      var = gsl_stats_variance_m(v over segment, m)
      madSeg = median( sorted |v_k - m| )              # MAD about the mean
      nfSeg  = madSeg / sqrt(m)                       # per-segment noise factor
      pooled running mean/var via tools::Statistics::collectedMean / collectedVariance   (Inferred: standard pooling)
      allMsd += madSeg
  allSd  = sqrt(pooledVar)
  allMsd = allMsd / countAccepted
  nf     = median(all nfSeg)
  return {mean: pooledMean, sd: allSd, msd: allMsd, nf: nf, n: L * countAccepted}
```

- Note: the n field is typed long long; `L*count` is Inferred from the cast pattern.
- Source of `s` (the derivative) and `pooledVar` details: Inferred.

### 2.2 AbstractPeakDetector<int,float>::calculateBaseline(Peak*) (MCP 0x6ed5ee; virtual, called as checkPeak via vtbl[4])

```
calculateBaseline(p):
  v[i] = EIC intensities, i = 0..n-1, scan index relative to firstScan
  minL = min over [0, maxScan-firstScan)   -> index minScanL
  minR = min over [maxScan-firstScan, lastScan-firstScan] -> index minScanR
      (both start at the apex value and index; strict <, so the first minimum wins)
  line through (minScanL, minL) and (minScanR, minR): slope = (minR-minL)/(minScanR-minScanL), icpt = minL - minScanL*slope
  abs = minScanL if minL < minR else minScanR (a tie takes the right point)
  repeat:
      for any i with (i*slope + icpt) - v[i] > 0.1: re-anchor the line through (abs, v[abs]) and (i, v[i])
  until no point violates
  residual r_i = v[i] - (slope*i + icpt) for all points, map<residual, index>
      (std::map insert: duplicate keys keep the first entry)
  take the first floor(map.size()/2)+1 entries in ascending residual order (size after duplicates collapse)
  linear fit gsl_fit_linear(x = index, y = v[index]) -> intercept, slope
  nu = sqrt(maxValue) * this->nf
  if nu == 0: warn "can not calculate noise units for s/n ratio. Rejecting this peak."; reject
  h  = maxValue - (intercept + slope*(maxScan - firstScan))
  sn = roundSig(h / nu, 8); setSignalToNoise(sn)
  accept iff h >= AN_MIN_PEAK_HEIGHT * nu
```

- Constant 0.1 in the re-anchor test: Observed.
- The baseline is a robust line through the lowest-residual half of the points: Inferred from the sort/fit pattern.

### 2.3 Noise factor

- AbstractPeakDetector ctor (MCP 0x6ed4f0): vptr set, nf = 1.0. Observed.
- Peaks created by MatchedPeakDetector get nf = 0.5 (Observed, Section 3.3).
- Peaks created by DefaultPeakDetector get nf = getNoiseFactor() (Observed). detectPeaks does NOT write analyzeBaseline's `nf` into this->nf: it has no setNoiseFactor call and no store to this+0x8 (Observed). The only setNoiseFactor caller is the GCMSAnalyzer<int,float> ctor (syms 0x5e88fc, call at 0x5e89bf); it passes a double from analyzer+0x70, filled by an indirect call at 0x5e8995: scan vtbl+0x50, `getNoiseFactor` → `DataProvider::getGlobalNoiseFactor` (Observed). On all 15 example runs it logs "Error determining nf. Using nf=1", so nf = 1 (Observed with the golden driver, `tools/golden/`).
- `GCMSDiskScan::determineNoiseFactor` (0x6e565a, objdump) fills the global statistics. For DataType 2 it calls `calculateGlobalBaselineStatistics(15, 50)`, which falls back to nf = 1. Otherwise it combines the per-scan `scan2stat` of the cache (file-formats.md §2.2): pooled mean and variance via `collectedMean`/`collectedVariance`, nf = median of the per-scan nf (field 4), n = sum of n; stored as `BaselineStatistics(mean, sqrt(var), 0, nf, n)` and written to the index. With zeroed per-scan statistics nf = 0 and every peak fails "can not calculate noise units" (Observed with the golden driver). With nf = 1 and mean = sd = 1 the old app also moved some peak edges against mean = sd = 0, so the pooled mean or sd feed the detector somewhere (Observed, not analysed).

## 3. Peak detection

### 3.1 DefaultPeakDetector<int,float>::detectPeaks (MCP 0x6ef020; syms 0x5ef020)

Setup (Observed):

```
if N < 5: return empty
BEGIN = AN_PEAK_THRESHOLD_BEGIN ; END = AN_PEAK_THRESHOLD_END
D  = SavitzkyGolayFilterFirstDerivative<float>(float(values))   # syms 0x64ab49, see below
S  = SavitzkyGolayFilter<int>(values)                         # syms 0x645577, see below
stat = analyzeBaseline(S, D, 15)                              # Section 2.1
```

Filters (Observed, 5 points each):

- `S[i] = (int) trunc((-3*(v[i-2]+v[i+2]) + 12*(v[i-1]+v[i+1]) + 17*v[i]) / 35.0)`; S[0], S[1] copy S[2] and S[N-2], S[N-1] copy S[N-3].
- `D[i] = (-2*v[i-2] - v[i-1] + v[i+1] + 2*v[i+2]) / 10.0f`; the first two and last two values are 0.
- `border_extension` is member this+0x18, set to false by the ctor (0x5eef5a); `setPeakBorderExtensionEnabled` (0x5ef9fc) has no direct caller (Observed).

Loop (reconstructed from objdump, syms 0x5ef1a7..0x5ef773; Observed control flow, names Inferred):

```
i = 2; tailing = false
while i < N - 2:
    beg = 0; mx = 0; en = 0
    if D[i] > BEGIN or tailing:
        beg = i; mx = i; zc = 0; i += 1
        while i + 1 < N:
            if D[i] >= 0:
                if D[i] == 0:
                    zc += 1
                    if zc > 2 and |S[i] - stat.mean| <= stat.sd: mx = 0; break
                mx = i; i += 1; continue
            else:
                if D[i+1] < 0: break
                i += 1; continue
        if mx == 0: tailing = false; continue
        while i + 1 < N:
            if D[i] > END and D[i+1] > END: en = i - 1; break
            i += 1
        if en == 0: tailing = false; continue
        if en - beg <= 3 or mx == beg or mx == en: tailing = false; continue
        p = createPeak(ion, beg, en, mx, c_values, D, border_extension, stat)   # Section 3.2
        if not checkPeak(p):                    # vtbl[4] = calculateBaseline (Inferred identity)
            tailing = false; i = beg
        else:
            if tailing: tailing = false
            t1 = c_values[lastScan] - stat.mean          # RAW intensities, absolute index
            t2 = c_values[maxScan]  - stat.mean
            if t1 > 0.2 * t2: tailing = true              # C = 0.2 at rodata 0xa96ca8 (value read from disassembly, Observed)
            peaks.push_back(p)
    i += 1
```

Notes:
- After an accepted peak the loop does `i -= 1` (0x5ef73d) before the final `i += 1`, so the next search starts at the scan where the falling flank ended; a tailing peak starts at the previous last scan + 1. Observed and confirmed by the golden test.
- Behavior when `i` equals the break index is Inferred.

### 3.2 DefaultPeakDetector<int,float>::createPeak (MCP 0x6efa1c; syms 0x5efa1c)

```
createPeak(ion, begin, end, max, values, D, extend, stat):
  if extend:
      while begin-1 valid and v[begin-1] < v[begin] (monotone toward baseline) and |v[begin] - stat.mean| >= 3*stat.sd: begin -= 1
      while end+1 valid  and v[end+1]  < v[end]   and |v[end]   - stat.mean| >= 3*stat.sd: end += 1
  apex = max; for k in begin..=end: if v[k] > v[apex]: apex = k   (strict >, starting from the flank search's max, so a tie keeps that scan; Observed, objdump 0x5efa1c)
  return Peak(ion, begin, end, apex, values[begin..end], nf = getNoiseFactor())
```

- The exact border conditions (direction check, and which side of stat is used) come from the decompile. Monotonicity check is Observed; exact comparison operator direction: Inferred.
- The "3 * sd" constant: Observed.

### 3.3 MatchedPeakDetector<int,float>::detectPeaks (MCP 0x6ee9ce; syms 0x5ee9ce)

Members (ctor 0x5ee926, Observed): `peak_width` (size_t) at this+0x10, `threshold` (int) at this+0x18. The base ctor sets nf = 1.0, which createPeak does not read.

```
density = 4 * peak_width                                   # shl $2 at 0x5eea0d (Observed)
if values.size == 0: return empty
f = tools::Filter::matched2ndDerivativeGaussianFilter<int>(values, density)
prev = f.at(0); i = 0
while i < f.size:
    begin = end = i
    if prev < thr and f[i] > thr:                          # upward crossing, both strict
        prev = f[i]; i += 1                                # prev is not updated in the inner loop
        while i < f.size:
            if prev > thr and f[i] <= thr: end = i; break  # first scan at or below thr is INCLUDED
            i += 1
    if begin != end and end - begin > 2:
        peaks.push_back(createPeak(ion, begin, end, *f))   # the FILTERED vector (Observed, rbx at 0x5eec6d)
    prev = f[i]; i += 1                                    # f[f.size] when the region stays open (out of bounds; the loop then ends)
```

All of this is Observed from objdump. A region still above the threshold at the trace end gets `end == begin` and is dropped. The next region can start at `end + 1` at the earliest, and only if `f[end] < thr`.

Filter (syms 0x64a5fd int, 0x64966d double; Observed): for x in [-d, d], u = 4x/d, kernel k(x) = (1.73/sqrt(d)) * (1 - u^2) * exp(-u^2/2), a negative 2nd-derivative Gaussian with sigma = d/4. Taps outside the signal count as 0; the int variant truncates the result (cvttsd2si). Both variants then overwrite the ends: `out[i] = out[d]` for i < d, then `out[i] = out[n-d-1]` for i >= n-d, or 0 everywhere if n <= d (0x64a901..0x64aa61 int, 0x64996d..0x649af7 double). So the zero taps never reach the result, and a filter without this step finds 333 extra peaks at the trace ends of 7235eg04 (matched golden comparison).

GUI (Observed): the only ctor caller is `CalibrationSetter::threadFinished` (0x4947f4), the "Compound Detection" tab of the *Calibration Setter* window. Its combo box "Peak detection method" offers "Default" and "Matched Filtering"; the choice is not saved in the settings file.

- "Default" (`currentText() == "Default"`, 0x494674): the five fields write the globals AN_PEAK_THRESHOLD_BEGIN (spin box +0x210), AN_PEAK_THRESHOLD_END (+0x228), AN_MIN_PEAK_HEIGHT (+0x240), AN_NUMBER_BINS (+0x258) and AN_DECONVOLUTION_WIDTH (+0x270, labelled "Peak sharpness proportionality:"). Then `new DefaultPeakDetector()`.
- Otherwise: `MatchedPeakDetector(peak_width = (int) +0x228, threshold = (int) +0x210)`. The other three fields are NOT written, so bins, deconvolution width and minimum peak height keep their current values.
- `on_combo_method_currentIndexChanged` (0x4938b0) relabels the first two fields. Default: "Peak threshold beg.:" = 10, "Peak threshold end:" = -10 (range -1e8..0). Matched: "Peak threshold:" = 3, "Peak width:" = 10 (range 1..1000). Values from rodata 0xa7e8b8..0xa7e908 (Observed).
- `CompoundDetector::run` (0x4952d8) then runs the same pipeline as the main window: `GCMSAnalyzer(scan, detector)`, `detectAllPeaks(true, 10.0f, 0.05f)`, `deconvoluteChromatogram()`, `calculateRI`, `Compound::toDisk` (.cmp).

Consequence (Observed): the deconvolution (Section 4) reads the peaks' own values through `Peak::getValue` / `getMaxValue`. With the matched detector these are the filtered values, so the apex histogram, the model peak, the fit and the raw integral all use the filtered traces.

### 3.4 MatchedPeakDetector<int,float>::createPeak (MCP 0x6eed8e; syms 0x5eed8e)

- apex = argmax over [begin, end] (strict >, starting at begin) on the supplied (filtered) values.
- Peak(ion, begin, end, apex, values[begin..=end], nf = 0.5). No `checkPeak`, so S/N keeps the PeakEssential default -1 (Observed in the golden dump).

### 3.5 Peak<int,float>::Peak(ion, begin, end, maxScan, values, nf) (MCP 0x6f2754; syms 0x5f2754)

- PeakEssential<float>(ion, begin, end, maxScan); stores values and nf.
- Throws if (end - begin + 1) != values.size() (Observed).
- If ion != 0: determinePreciseMaximizationTime(); determinePeakQuality(); calculatePeakSharpness().

### 3.6 determinePreciseMaximizationTime (MCP 0x6f3782; syms 0x5f3782)

```
points: scans m-1, m, m+1 around maxScan (m = apex)   # Observed
X_i = AN_NUMBER_BINS * scan_abs_i ; Y_i = v_i
(c0, c1, c2) = tools::Statistics::polynomialRegression(X, Y, 2)   # Observed
t = round(-c1 / (2*c2)) / AN_NUMBER_BINS                       # rounding: Observed as (x + 0.5) truncated; sign conventions Inferred
if |t - maxScan| < 1.0: setInterpolatedMaxScan(t)
else: warn "Timeshift is greater or equal 1"; no interpolated max
bounds failure: "Strange peak ... Can not calculate precise maximization time!"
```

`polynomialRegression(X, Y, 2)` (0x6c3292) builds the Vandermonde matrix `[1, x, x²]` and calls `gsl_multifit_linear` (GSL 2.x: column balancing, SVD, truncation at `DBL_EPSILON`). The rounding is `(int)(t + 0.5)` (`cvttsd2si` at 0x5f3b0b). With integer intensities the exact vertex often lies on a half step. Which way the old code rounds such a tie then depends on the fit's rounding error in absolute bin coordinates. On the 15 example runs it rounds 217 of 267 ties down, and three equal values (flat top) give an arbitrary vertex near the apex. Observed with the golden driver (issue 46).

### 3.7 determinePeakQuality (MCP 0x6f3d28; syms 0x5f3d28)

```
norm = 0; pen = 0
for i in 0 .. n-2:
    s = v[i+1] - v[i]
    if firstScan + i < maxScan:           # left of apex
        if s >= 0: norm += s  else: pen += |s|
    else:                                 # at or right of apex
        if s <= 0: norm += |s| else: pen += s
quality = 100 if norm == 0 else 100 * pen / norm
if quality is NaN: print "Bad quality..." and do not set
```

- Observed.

### 3.8 calculatePeakSharpness (MCP 0x6f3f98; syms 0x5f3f98)

```
nu = sqrt(maxValue) * nf
sharpL = max over i < apex with v[i] >= 0 of (maxValue - v[i]) / ((apex - i) * nu)
sharpR = max over i > apex with v[i] >= 0 of (maxValue - v[i]) / ((i - apex) * nu)   # symmetric; exact form Inferred
sharpness = max(0, (sharpL + sharpR) / 2), rounded roundSig 8
```

- Observed, both sides with the same `v[i] >= 0` test. Quirk (Observed): sharpL and sharpR start at -1.0, so a side without candidates adds -1 before the clamp; NaN/inf results are an error path.

## 4. Deconvolution / component detection

All addresses below are syms addresses (MCP = +0x100000). Ghidra cuts most of these functions at `std::endl` and `operator new` (it marks both no-return). The control flow after those calls was read from objdump.

### 4.0 Pipeline and entry points (Observed)

- GUI path: `metabolitedetector::ComponentDetector::detectComponents()` (0x4b3a0c) calls `GCMSAnalyzer::detectAllPeaks(reorder=true, scan_dist=10.0f, ion_dist=0.05f)` (call at 0x4b3a87; floats at rodata 0xa81010 = 10.0 and 0xa8100c = 0.05). Then it calls `deconvoluteChromatogram()` (0x5e907e).
- RI calibration wizard: `RICalibrationWizard::CompoundDetector::detect()` calls `detectAllPeaks(reorder=false, ...)` (0x3f4710), then `deconvoluteChromatogram()`.
- Batch quantification: `ResultSeries::scanIons` (0x638589..0x6385d7) builds an analyzer with a fresh `DefaultPeakDetector`, calls `setMinimumPeakNumberRequired(1)`, then `deconvoluteChromatogram(peaks, scan_count, start_scan, 2 * AN_DECONVOLUTION_WIDTH)` on a scan window. Details of the whole batch rescan: [batch-quantification.md](batch-quantification.md).
- `MIN_PEAK_NUMBER` is field +0x0 of GCMSAnalyzer. The settings dialog sets it from `AN_MIN_PEAK_NUMBER` (call at 0x3c1177). The GCMSAnalyzer ctor sets it to `AN_MIN_PEAK_NUMBER` (Observed).
- `deconvoluteChromatogram()` (0x5e907e) = `deconvoluteChromatogram(this->peaks, scan->vtbl[3]() /*scan count*/, 0, AN_DECONVOLUTION_WIDTH)`. Observed.

### 4.1 detectAllPeaks(reorder, scan_dist, ion_dist) (0x5e8ac0)

```
clear peaks, bins, components, compounds
for i in 0 .. scan->vtbl[4]() - 1:              # number of ion traces
    chrom = scan->vtbl[5](i)                     # EIC of ion i
    peaks.splice(end, detector->detectPeaks(chrom))   # Section 3
    progress every 100 ions; abort on cancel -> return false
print "detected N peaks..."
if cancelled: return false
if reorder and scan->getDataType() != 2: reorderPeaks(scan_dist, ion_dist)
return true
```

Observed (tail at 0x5e8e99..0x5e8f15). Deconvolution passes 1.0 as the Compound "resolution" for type 2, else 0.0 (Section 4.6).

Data type, mass traces and centroiding (Observed):

- Data type 2 = nominal masses. The CDF import (`CdfImporter::scanAdd`, 0x6cd038) always writes type 2. Only `FastFlightImporter::importData` (0x6c4cde) writes type 1 or 0 and keeps exact masses. So reorderPeaks (4.2) runs only for FastFlight runs.
- `FastFlightImporter::importData(path, drop_first, b)` writes DataType `b ? 1 : 0`; `b` is the user option at `FastFlightBatchImporter`+0x60 (Observed). Type 1 caches are compressed profiles with per-scan baseline statistics (file-formats.md §2.2). `GCMSDiskScan::getChromatogram` reads type 0 through `getProfileChromatogram` (not ported; on our caches it throws "start must be lower than stop"), types 1 and 2 through the ion list below.
- `scanAdd` per scan:
  - Masses become float (`cvtsd2ss`) and are multiplied by the mass scale factor when it is not 1.
  - Intensities are scaled, then rounded with `(int)(x + 0.5)`.
  - Then `Spectrum<int,float>::getCentroidSpectrum(0)` (0x600ffc, objdump) runs on the scan.
- getCentroidSpectrum, from mass index i on:
  - A part continues until the m/z gap to the next point is > 0.1 (rodata 0xa96cb0).
  - Or until a valley: first two falls in a row (`v[i+1] < v[i]` and `v[i+2] < v[i+1]`), then two rises in a row. Past the end, the missing neighbour repeats the last value read.
  - The part ends at that point; the next part starts after it.
  - max = the largest value in the part. A part of one point is kept as it is.
  - Otherwise: `intensity = num = denum = 0`. Every point with `v >= SF_FWHM (0.5) * max` adds v to intensity and denum, and `(float)v * m` to num. Then `centroid = num / denum`, or -1 if both are 0.
  - With type 1 the intensity is max instead (the importer passes 0).
  - Parts with a negative centroid are dropped. Value = `(int)(intensity + 0.5)` (`is_fl_value` is 0). Mass = `(float)centroid` (`is_fl_mass` is 1).
- The masses in the old disk-scan cache are these centroids rounded to nominal masses, `floor(m + 0.5)`. Checked against the cache of the synthetic high-resolution run: all 4400 scans are equal. Where exactly the rounding happens is not analysed.
- On unit-mass data every gap is 1, so each point is its own part; centroiding then only rounds intensities and drops negative masses.
- `GCMSDiskScan::getChromatogram` reads each trace with `getSignal(ion)` (0x5e0aec): the point of exactly that float mass.
  - A missing mass reads `(int)` of the scan's first baseline statistic (0 for CDF), or 0 when baseline correction is on.
  - A found point is baseline-corrected (vtbl+0x48) when that is on.

Port: `IonTraces::new` centroids like this (`centroided`; skipped when no scan has two points within 0.1, then fractional intensities and m/z in [-0.5, 0) stay as read) and sums the centroids per nominal mass. `IonTraces::exact` makes one trace per distinct m/z, with 0 where a scan lacks that mass; `reorder_peaks` is 4.2. `Chromatogram::exact_masses` marks runs with exact masses (FastFlight, caches with DataType 0 or 1); `IonTraces::new` then builds exact traces and `IonTraces::detect_all_peaks` runs `reorder_peaks(10, 0.05)`, like `detectAllPeaks(true, 10, 0.05)`. Components keep the exact masses in `Component::exact_spectrum` (the old `map<float,int>`). Golden check: `old-md-3.3-exact-masses.tsv`, a type 1 cache of the synthetic high-resolution run with per-scan statistics (0, 0, 0, nf = 1, 100).

### 4.2 reorderPeaks(scan_dist, ion_dist) (0x5eac2c) - align apexes of split masses

Used for non-type-2 (FastFlight, exact-mass) data, where one nominal fragment can be split over neighbouring m/z traces. With ion_dist = 0.05 it never merges unit-mass traces. Observed (Ghidra decompile, structure confirmed):

```
peaks.sort(Peak_sort(ION))                      # by ion asc, then firstScan asc
all_ions = set(scan->getAllIons())
NB = AN_NUMBER_BINS
bin_dist = (uint) roundf(NB * scan_dist)
bins = vector<list<IonProfile*>>(NB*(scan_count+1) + 1)
for p in peaks:
    b = (uint)(p.interpolatedMaxScan * NB + 0.5)       # no first_scan offset here
    if b >= bins.size: throw
    if bins[b] empty: insertIonProfile(b, bins, bin_dist, p); continue
    winner = null; min_dist = INT_MAX
    for prof in bins[b]:
        nxt = successor of prof.last_ion in all_ions   # throws if last_ion not in set; warns if none
        if p.ion == nxt and |prof.last_ion - nxt| <= ion_dist:
            d = |p.maxScan - prof.max_scan|
            if d < min_dist: min_dist = d; winner = prof
    if winner == null: insertIonProfile(b, bins, bin_dist, p)
    elif p.firstScan < winner.max_scan < p.lastScan: winner.addPeak(p)
    # else: p is not put in any profile
for each profile (each once, `processed` flag):
    if profile has >= 2 peaks:
        s = Peak(profile.peaks)                         # summed peak, 4.3
        for q in profile.peaks with q.firstScan <= s.maxScan <= q.lastScan:
            q.setMaxScan(s.maxScan); q.setInterpolatedMaxScan(s.interpolatedMaxScan)
delete all profiles
```

- `insertIonProfile(bin, bins, dist, p)` (0x5eb972, objdump): `prof = new IonProfile(p)`. Then `lo = (int)(bin - dist) < 0 ? 0 : bin - dist` and `hi = min(bin + dist, bins.size - 1)`. It pushes `prof` into `bins[k]` for k in [lo, hi). The upper bound is exclusive, so bin `hi` is not filled (Observed: `jae` at 0x5eba59). It appends (`list::push_back`, 0x5eba7b), and reorderPeaks walks `bins[b]` from `begin()` with `operator++` (0x5eaf62..0x5eb1b3). With the strict `d < min_dist`, the profile created first wins a tie (Observed).
- `IonProfile(p)` (0x649594): processed = false; addPeak(p). `IonProfile::addPeak(p)` (0x648f90): push p; last_ion = p.ion; max_scan = Peak(peaks).maxScan (the apex of the summed profile). Observed.

### 4.3 Peak(list<Peak*>) - summed peak (0x5f28f8, Observed)

```
begin = min firstScan; end = max lastScan; sn = mean of member S/N
values[i] = sum over members covering scan i of member.getValue(i)   # int sum
maxScan = first argmax(values); throws if the list is empty
determinePreciseMaximizationTime(); determinePeakQuality(); calculatePeakSharpness()   # Section 3.6-3.8
```

`Peak::getValue(scan)` throws if scan is outside [firstScan, lastScan] (0x5f2efa). `Peak::getMaxValue()` = values[maxScan - firstScan], the raw apex intensity (0x5f2ea6).

### 4.4 initializeBins(scan_count, first_scan) (0x5ec31c, Observed)

```
NB = AN_NUMBER_BINS
bins = vector<Bin>(NB * (scan_count + 1) + 1); bins[i].ID = i
for p in peaks:
    b = (uint)(NB * (p.interpolatedMaxScan - first_scan) + 0.5)
    if b >= bins.size: print "Invalid interpolated max scan ..." and skip p
    else: bins[b].addPeak(p)
Bin::addPeak(p) (0x649e58): bin.peaks.push_back(p); bin.overall_sharpness += p.getMaxValue()
```

So the field called `overall_sharpness` is the SUM OF RAW APEX INTENSITIES in the bin (a binned "apex TIC"), not a sharpness. Observed. AN_NUMBER_BINS is the number of sub-scan bins per scan. It also sets the time quantum of the precise apex in 3.6.

### 4.5 detectComponents(width) (0x5ebac4, objdump; Observed)

```
components.clear()
W = (size_t)(AN_NUMBER_BINS * width)                       # width = AN_DECONVOLUTION_WIDTH (scans)
v[i] = bins[i].overall_sharpness
f = tools::Filter::matched2ndDerivativeGaussianFilter<double>(v, W)   # kernel: Section 3.3
T = 5.0                                                     # rodata 0xa96c60
prev = (int) f[0]; i = 0
while i < f.size:
    begin = end = i; t_peaks = []
    if prev < T and (int) f[i] > T:                         # upward crossing (on int-truncated values)
        t_peaks += bins[i].peaks; hold = (int) f[i]          # hold > T, never updated
        loop:
            i += 1; if i >= f.size: break
            t_peaks += bins[i].peaks; end = i
            if hold > T and (int) f[i] <= T: break           # first bin at/below T is INCLUDED
    if begin != end:
        n = begin                                            # extend left down the flank
        while n > 1 and f[n] > f[n-1]: t_peaks += bins[n].peaks; begin = n; n -= 1
        n = end                                              # extend right down the flank
        while n < f.size-1 and f[n] > f[n+1]: t_peaks += bins[n].peaks; end = n; n += 1
        components.push(DeconvolutedCompound(begin, end, t_peaks))
    prev = (int) f[i]; i += 1
```

Quirks that a faithful port must keep or consciously drop (all Observed):
- Both flank extensions start at n = begin / n = end. They therefore insert `bins[begin].peaks` (if f[begin] > f[begin-1]) and `bins[end].peaks` (if f[end] > f[end+1]) a SECOND time. `t_peaks` then holds duplicate Peak pointers. This changes the peak count used against MIN_PEAK_NUMBER, the 25 % model-peak quota, the ModelPeak sum (a duplicated model peak counts twice) and the stored compound peak list (4.6).
- The flank extension stops one bin before the local minimum of f. The stop bin itself is not added.
- A region of a single bin (begin == end) is dropped.
- Comparisons use `(int)` truncated filter values against T = 5.0.
- `DeconvolutedCompound(start_bin, stop_bin, peaks)` (0x5ed4a2) only stores the three values.

### 4.6 deconvoluteChromatogram(peaks, scan_count, start_scan, width) (0x5e90de, Ghidra + objdump; Observed)

```
this.peaks = peaks; compounds.clear()
initializeBins(scan_count, start_scan); detectComponents(width)
for comp in components:
    t = copy(comp.peaks)
    if t.size < MIN_PEAK_NUMBER: print "Rejecting component because only N peaks were detected."; continue
    sort t by sharpness DESC                  # Peak_sort(SHARPNESS): p2.sharp < p1.sharp (0x6a8252)
    quota = (uint)(t.size * 0.25); if quota == 0: quota = 1
    m = []
    for p in t:                               # pick model peaks
        if p.quality <= 10.0 and p.ion not in any [lo, hi] of QUANT_EXCL_IONS:
            if quota == 0: break
            m.push(p); quota -= 1
    if m empty:
        print "WARNING: Using peak of bad quality for model peak construction..."
        s0 = t[0].sharpness * 0.75            # rodata 0xa96c90
        for p in t: if s0 > p.sharpness: break; m.push(p)      # prefix of the sharpness-sorted list
    model = ModelPeak(m)                      # 4.7
    ions = []; ints = []
    for p in t (erase-while-iterating):
        try: c = calculatePeakFunction(p, model)              # 4.8
        except GCMSScanException: print "Error during model peak function calculation: ..."; erase p; continue
        if c empty: erase p; continue
        h = (int)(c[2] * model.getMaxValue())                 # modelled apex height of this ion
        if h < 1: erase p; continue
        ions.push(p.ion); ints.push(h)
        calculatePeakIntegral(model, c[2], p)                 # 4.9, stores p.integral
        p.setModelParameter(c[2])
    if ions.size < MIN_PEAK_NUMBER: print "Less than N peaks for this compound. Skipping this compound."; continue
    if ions.size != ints.size: print "Sizes do not match."; continue
    mapping = map<float,int>(); bp = 0
    for k: mapping.insert(ions[k] -> ints[k]); bp = max(bp, ints[k])   # duplicate ion: FIRST (= sharpest) wins
    minh = AN_MIN_HEIGHT_BP * bp                            # 0x5e9bf2 (double, value from settings)
    spectrum = [(ion, h) in mapping (ascending ion) with h >= minh]    # others erased from mapping
    if spectrum empty: print "No values(sorted_ions)... During deconvoluition..."
    res = (scan->getDataType() == 2) ? 1.0 : 0.0
    rt  = scan->getRetentionTime((double) model.getMaxScan())
    C = new Compound(rt, &ints_sorted, &ions_sorted, (float) res); C.setScan(model.getMaxScan())
    C.setPeaks([copy of PeakEssential(p) for p in t if p.ion in mapping])   # duplicates kept
    C.setModelPeak(Peak<double,short>(ion=1, model.begin, model.end, model.maxScan, model.values, nf=0.5))
    compounds.push(C)
    if C.values empty: print "No values... During deconvoluition..."
scan->getRICalculator().calculateRI(compounds)          # Section 6
```

- The QUANT_EXCL_IONS test is inclusive on both ends (`lo <= ion <= hi`). Observed.
- AN_MIN_HEIGHT_BP is a fraction of the base peak of the deconvoluted spectrum. Its default value is 0 (`.bss`, Section 1).
- The compound spectrum is the set of modelled apex heights `c2 * max(model)`, rounded toward zero. It is NOT the raw apex intensities. Observed.

### 4.7 ModelPeak(vector<Peak*> m) (0x647ce4, Observed)

```
if m empty: throw
begin = min firstScan; end = max lastScan
values[i - begin] = sum over p in m with p.first <= i <= p.last and p.size > 0 of (double) p.getValue(i)
maxScan = first argmax of values (strict >), scanning begin..end
getValue(scan) = values[scan - begin], prints "scan, begin, end" to cerr if outside (0x64852c)
getMaxValue() = getValue(maxScan)   (0x648176)
```

The model peak is the plain SUM of raw intensities of the selected model peaks. It is not normalised.

### 4.8 calculatePeakFunction(peak, model) (0x5ea4f0, Observed)

Ordinary least-squares fit of each ion trace against the model shape over the overlap of the two scan ranges:

```
lo = max(model.begin, peak.firstScan); hi = min(model.end, peak.lastScan)
n = hi - lo + 1; if n < 3: return []
X row k = [1, (double) scan, model.getValue(scan)], Y_k = peak.getValue(scan), scan = lo + k
gsl_multifit_linear(X, Y) -> c = (c0, c1, c2)
return [roundSig(c0, 8), roundSig(c1, 8), roundSig(c2, 8)]
```

So intensity_ion(scan) ~ c0 + c1*scan + c2*model(scan). c0 + c1*scan is a local linear baseline. c2 is the share of the model peak in this ion trace. `roundSig` is `PeakEssential::roundSignificantFigures(x, 8)`.

### 4.9 calculatePeakIntegral(model, c2, peak) (0x5ea99e, Observed)

```
lo = max(model.begin, peak.firstScan); hi = min(model.end, peak.lastScan)
integral = sum over scan in [lo, hi] of (AN_MODELED_VALUES_FOR_INTEGRATION ? c2 * model.getValue(scan) : peak.getValue(scan))
peak.setIntegral(integral)
```

Plain sum over scans (no trapezoid, no baseline subtraction, no time scaling). Observed.

### 4.10 Settings summary (deconvolution)

| Setting | Use | Evidence |
|---|---|---|
| AN_DECONVOLUTION_WIDTH | width W (in scans) of the matched 2nd-derivative Gaussian filter on the binned apex-intensity trace; filter length = AN_NUMBER_BINS * W bins | 0x5e907e, 0x5ebb30 |
| AN_NUMBER_BINS | sub-scan bins per scan (bin = round(NB * (interpMaxScan - first_scan))) and precise-apex quantum | 0x5ec31c, 0x5ebb07 |
| AN_MIN_PEAK_NUMBER | minimum ion peaks per component, before and after fitting | 0x3c1177, 0x5e9248, 0x5e9648 |
| AN_MIN_HEIGHT_BP | relative cutoff of modelled ion heights vs. compound base peak | 0x5e9bf2 |
| AN_MODELED_VALUES_FOR_INTEGRATION | integrate the model fit instead of raw data | 0x5ea99e |
| QUANT_EXCL_IONS | ions excluded from model-peak selection only | 0x5e946a |

### 4.11 Unknowns (deconvolution)

- (resolved) `tools::Filter::matched2ndDerivativeGaussianFilter<double>` edge handling: Section 3.3.
- Meaning of `getDataType() == 2`, and how `Compound(rt, ints, ions, res)` uses `res`.
- `scan->vtbl[3]` (scan count) and `vtbl[4]/[5]` (ion count / EIC) identities are Inferred from usage.
- Whether the duplicated bin peaks from 4.5 are intended. A port should reproduce them first to match golden output.

## 4b. Centroided and nominal spectra for the NetCDF export (Observed; ported in `crates/md-core/src/centroid.rs`)

The old `NetCDFExporter` writes each scan as profile, centroid or nominal mass (`_profile/_centroid/_nominal.cdf`) with
`NetCDFExport::addScan(ions, intensities, tic, rt_ms / 1000)`. The TIC is the scan's own, not the sum of what is written.
The ports below are checked bit for bit against the old binary's own functions: `tools/spectrum/run.sh` calls
`SpectrumFunctions::getCentroidSpectrum<int,float>` (0x5dfe78) and `Spectrum<int,float>::getBinnedSpectrum` (0x602a6a) on
4000 random spectra of 3 to 600 points (3180 centroids); every mass, area and point count equals the port's
(`crates/md-core/tests/centroid.rs`, `matches_old_program`).

### 4b.1 getCentroidSpectrum<int,float>(spectrum, begin, end, maxdev) -> detectPeaks<int>

Thresholds are `SF_CENTROID_THRESHOLD_BEGIN/END/MAX_BASE_DEV` (10, -5, 30). Values are the ints of the spectrum.

```
d[i] = trunc((-2 v[i-2] - v[i-1] + v[i+1] + 2 v[i+2]) / 10.0) for 2 <= i < n-2, else 0   (SavitzkyGolayFilterFirstDerivative<int>)
n <= 4: no peaks
(mean, sd) = analyzeBaseline<int>(v, 25)
i = 2; while i < n-2:
  if d[i] > begin:
    start = top = i; i++
    while i < n-2: if d[i] >= 0 { top = i; i++ } else if d[i+1] < 0 { i++; break } else { i++ }
    end = 0
    while i+1 < n:                         # the same i: the next peak is searched after this position
      if d[i] > endthr && d[i+1] > endthr && (sd == 0 || sd*maxdev > |v[i-1]-mean|) { end = i-1; break }
      i++
    keep if end != 0 && end-start > 2 && top != start && top != end
  i++
```

A kept peak covers `v[start..=end]` (the extension by the peak width exists in the old code but the centroid caller leaves it off):
- `createSpectrumPeak`: highest point `max` = first strict maximum from `top`.
- area (`determineIntegral`, called before the baseline is set): the plain sum of the points, truncated to int. `determinePeakBaseLine` runs afterwards and does not change it.
- mass (`determineCentroid`): over the points with `v >= 0.5 * v[max]`: `x += (float)mz * (float)v` (float product, summed in double), `w += v`, count += 1; mass = `(float)(x / w)`. The count is the "FWHM".

`analyzeBaseline<int>(v, 25)`: segments of 25 points, `i = 0; i < n - 25; i += 25`. A segment counts if its values cross its own mean at least `25 / 2` times (a crossing: `v[k] >= m && v[k+1] < m` or `v[k] <= m && m < v[k+1]`); then its mean and `gsl_stats_variance_m` are pooled (`tools::Statistics::collectedMean/collectedVariance`, sample variance of all accepted points). A spectrum shorter than 25 points makes the old code use segments of one point, whose variance is 0/0: sd is NaN and no peak end is ever accepted. A spectrum with no accepted segment has mean 0 and sd 0.

### 4b.2 Spectrum<int,float>::getBinnedSpectrum(values, ions, width) (nominal mass)

```
half = width / 2.0 (double);  bin(m) = (float)(int)((m + half) / width) * width
cur = bin(ions[0]); best = values[0]
for i in 0..n:
  if ions[i] >= cur + half: push (cur, best); best = values[i]; cur = bin(ions[i])
  else if values[i] > best: best = values[i]
  if i == n-1: push (cur, best)
```
The export uses width 1. An empty spectrum gives an empty spectrum ("Creating empty spectrum").

### 4b.3 Not ported

- `Spectrum::getCentroidSpectrum(type)` (0x600ffc) chooses between these two by `CentroidIntensityType` and data type: data that is already centroid (DataType 2) is copied. The rewrite has no data type, so a centroid export of such a run centroids it again.
- The old export also takes `reorderPeaks` for non-CDF runs (§4.2, `IonTraces::exact`); that is a separate analysis step, not part of the export.

## 5. Library search and scoring

Template instance `<int,float>`. vtable `LibrarySearch<int,float>` at 0xcb7258 (Observed):

| Slot | Function |
|---|---|
| 0 | getLibraryHits 0x62657a |
| 2 | calculatePureSpecSimScore 0x621fa6 |
| 3 | calculateImpureSpecSimScore 0x6226d8 |
| 4, 5 | ...Wave variants 0x622b92 / 0x623e1a |
| 7 | filterMass 0x626236 |
| 11 | mapCompoundsToLibrary 0x629526 |

`Compound<int,float>` vtable (0xca2cc0, Observed):

| Slot | Function |
|---|---|
| 2 | getRetentionTime |
| 3 | getRetentionIndex |
| 6 | getCentroidCompound |
| 7 | getBinnedCentroidCompound(w) |
| 8 | getPeakEssentials |
| 16 | autoQuantificationIonAssignment(excl, min_dist, n, q) |
| 17 | autoQuantificationIonAssignment(n) |

### 5.1 Intensity weighting: Library::scaleCompound(Spectrum<int,float>) (0x621960, objdump; Observed)

```
for each peak i: scaled[i] = (float)( pow(ion_i, LS_SCALE_FACTOR_ION) * pow((double)intensity_i, LS_SCALE_FACTOR_VALUE) )
return new Spectrum<float,float>(ions, scaled) with the same bin width, RI and RT
```

This is the Stein-Scott style mass/intensity weighting m^a * I^b. LS_SCALE_FACTOR_ION is at 0xcc2138 and LS_SCALE_FACTOR_VALUE at 0xcc2140. No normalisation happens here; the cosine in 5.3/5.4 is scale-invariant.

### 5.2 Mass filter: Library::filterMass(ion) (0x626236, Observed)

Returns true (ion excluded) if `lo <= ion <= hi` for any pair in `LS_MASS_FILTER` (a set<pair<float,float>>, so it is sorted). It returns false early at the first pair with `ion < lo`. Excluded ions are skipped on the library side and are not counted as unmatched observed ions.

### 5.3 Pure score: calculatePureSpecSimScore(obs, lib) (0x621fa6, Ghidra; Observed)

Both inputs are scaled spectra (5.1). Ions are assumed ascending (Inferred: the merge logic needs it). This is a weighted cosine over the union of peaks. Unmatched peaks on either side are down-weighted by `NM = LS_SCALE_FACTOR_NO_MATCH`.

```
num = d1 = d2 = 0; io = 0; it = 0           # io: obs match cursor, it: obs "consumed" cursor
for k in 0 .. lib.n-1:
    L = lib.ion[k]
    o = io < obs.n ? obs.ion[io] : -1.0
    if filterMass(L): continue
    while io < obs.n - 1 and o < L: io += 1; o = obs.ion[io]
    dist = DBL_MAX
    if o != L:
        dist = |o - L|
        if io > 0 and |obs.ion[io-1] - L| <= dist: io -= 1; o = obs.ion[io]; dist = |o - L|
    while it < io and it < obs.n:            # observed peaks passed over = unmatched
        if not filterMass(obs.ion[it]): d2 += (NM * obs.val[it])^2
        it += 1
    vl = lib.val[k]; vo = 0
    if o == L or dist <= 0.01:
        vo = obs.val[io]; it += 1 (if it < obs.n); io += 1
    else:
        vl = NM * vl
    num += vl*vo; d1 += vl^2; d2 += vo^2
for remaining it < obs.n: if not filterMass(obs.ion[it]): d2 += (NM * obs.val[it])^2
return (d1 > 0 and d2 > 0) ? num / (sqrt(d1) * sqrt(d2)) : 0
```

- The match tolerance is a hard-coded 0.01 m/z (Observed literal), in addition to exact equality.
- Quirk (Inferred from the control flow): if step-back (`io -= 1`) lands on an observed peak that was already matched, that peak matches twice. `it` is then one ahead, and the next observed peak is never added to `d2`. This only occurs with closely spaced high-resolution ions.

### 5.4 Impure score: calculateImpureSpecSimScore(obs, lib) (0x6226d8, Ghidra; Observed)

Same matching as 5.3 without the `it` bookkeeping. Observed peaks that match no library ion are ignored, which gives a "reverse search" cosine over library ions only:

```
for each non-filtered lib ion: match as in 5.3 (tolerance 0.01)
   matched: vo = obs.val[io]; io += 1      unmatched: vl *= NM; vo = 0
   num += vl*vo; d1 += vl^2; d2 += vo^2
return (d1 > 0 and d2 > 0) ? num / (sqrt(d1)*sqrt(d2)) : 0
```

### 5.5 RI score: Library::calculateRISimScore(ri1, ri2, sigma) (0x622ae0, Observed)

```
if |ri1 - ri2| > 3*sigma: return 0
return exp( -(ri1 - ri2)^2 / (2*sigma^2) )          # sigma = LS_RI_DIFF
```

### 5.6 LibrarySearch::getLibraryHits(comp, bin_width) (0x62657a, Ghidra + objdump; Observed unless marked)

```
PF = LS_PURE_FACTOR (0xcc2120); IF = LS_IMPURE_FACTOR (0xcc2128); sigma = LS_RI_DIFF (0xcc2118)
has_ri = comp.RI >= 0
ri_c = has_ri ? comp.RI : comp.correctedRT / 1000
b_comp = comp
if comp.binWidth < bin_width:
    print "INFORMATION: Binning compound to bin width <w>"
    b_comp = comp.getBinnedCentroidCompound(bin_width)                 # vtbl[7], 0x6266f8
for L in lib (std::set<LibraryCompound*>, i.e. pointer order):
    if L.RI == -1.0 or not has_ri: ri_l = L.correctedRT/1000; ri_c = comp.correctedRT/1000 (if has_ri)
    else: ri_l = L.RI; ri_c = comp.RI
    ri = calculateRISimScore(ri_l, ri_c, sigma)
    if LS_USE_RI and ri <= 0: continue
    bring b_comp and L to a common bin width (rebinning the finer one with getBinnedCentroidCompound)
    spec_lib = (L unchanged) ? getPreScaledSpectrum(L) : scaleCompound(L_rebinned)   # pre-scaled cache, 0x62d428
    spec_comp = scaleCompound(b_comp) (cached while the bin width does not change)
    impure = vtbl[3](spec_comp, spec_lib); pure = vtbl[2](spec_comp, spec_lib)
    s = PF*pure + IF*impure
    o = 1.0
    if LS_USE_SPEC: o = s*s
    if LS_USE_RI:   o = o*ri
    if LS_USE_SPEC and LS_USE_RI: o = o^(1/3)      # = s^(2/3) * ri^(1/3)
    elif LS_USE_SPEC: o = sqrt(o)                   # = |s|
    hits.push(LibraryHit(L, spec_sim=s, ri_sim=ri, overall=o))
std::sort(hits) by overall ascending (LibraryHit::operator< 0x62dcbe); std::reverse -> descending
```

- If neither score is enabled, every hit has overall = 1.0 (0xa96c48). If only RI is enabled, overall = ri.
- When either RI is missing, the "RI" comparison uses corrected RT / 1000 with the same sigma (Observed). The time unit is Inferred to be ms, so the comparison is in seconds.
- The bin-width branch at 0x626800..0x626a80 has three cases: lib coarser than compound, lib finer than compound, and equal. Their outcome is summarised above (Inferred intent). The exact branch conditions are in the Ghidra output and were not re-verified in objdump.
- The `LibrarySearch(vector, bool)` ctor builds the pre-scaled cache (createScaledLib 0x62d10c → add2ScaledLib for every entry). The role of LS_USE_SCALED_LIB in that choice is Unknown.
- `AMSLibrary<int,float>::getLibraryHits` (0x62e8f4) is a variant: no rebinning or scaling, filter `ri >= 0` (always true), `overall = s * ri` with no cube root (Observed).

### 5.7 mapCompoundsToLibrary(comps, cutoff, alternatives*) (0x629526, Ghidra; Observed)

This is a greedy one-to-one assignment between detected compounds and library entries.

```
alternatives.clear(); unassigned = set(comps); all_hits = map<LibraryCompound*, list<LibraryHit>>
for comp in comps (pointer order):
    hits = getLibraryHits(comp, (float) LS_RESOLUTION)           # LS_RESOLUTION 0xcc2150 = bin width
    for h in hits while h.overall >= cutoff:                     # hits are sorted, so this stops at the first miss
        alternatives[comp].push(h)
        all_hits[h.libCompound].push(LibraryHit(comp, h.spec_sim, h.ri_sim, h.overall))
for each L: sort all_hits[L] by overall desc
while all_hits not empty:
    best = -1.0; best_ref = null
    for L in all_hits (pointer order):
        any = false
        for h in all_hits[L]:
            if h.comp in unassigned:
                any = true
                if h.overall > best: best = h.overall; best_ref = L; best_hit = h
        if not any: drop L
    if best_ref: result[best_ref] = best_hit; unassigned.remove(best_hit.comp); all_hits.erase(best_ref)
return result        # map<LibraryCompound*, LibraryHit(detected compound)>
```

Ties are broken by pointer order (first strictly greater wins), so tie results are not reproducible across runs (Inferred).

### 5.8 GUI flow: ComponentDetector::librarySearch (0x4b3bc2, Observed)

- It calls `mapCompoundsToLibrary(comps, cutoff = 0.0, &all_hits)`; the `pxor xmm0` at 0x4b3d2f shows the cutoff is zero.
- `compounds[comp]` = all hits for comp (descending). `mapping[comp]` = the assigned library hit, or `LibraryHit(null)` if there is none.
- `LS_CUTOFF` is not applied during the search. It is compared against `getOverallScore()` later:
  - DataTableModel::getData (0x4c88db)
  - CompoundMapper::findCompound (0x4b6e1e..0x4b6e32, hit accepted if overall >= cutoff; Inferred from the comisd/jb)
  - BatchQuantificator::startAnalysis/createLibrary

### 5.9 Unknowns (library search)

- Defaults of all LS_* settings.
- The exact `Spectrum::getBinnedSpectrum` / `getCentroidCompound` algorithms (0x6019d6, 0x60e212). These were not decompiled.
- The Wave score variants (vtbl 4/5) and `calculateScore_frag` (vtbl 6): when they are used.
- `LS_USE_SCALED_LIB` semantics. (`LS_IDENTIFICATION_ALGORITHM` is always 0: file-formats.md §5.)
- SQLLibrary::getLibraryHits (0x6e306c) reads LS_CUTOFF directly. Not read.
- Whether spectra are guaranteed sorted by ion before 5.3/5.4. They are built from std::map in 4.6, so deconvoluted spectra are sorted (Observed); library spectra are Unknown.

## 6. RI calibration and RT shift correction

### 6.1 RICalculator data (Observed)

- `rt2ri` and `ri2rt` are `std::map<double,double>`. `drt` is a global retention-time shift. `source` is a string (the ctor argument).
- `addCalibrationPair(rt, ri)` (0x6d21f0) inserts rt→ri and ri→rt. std::map::insert keeps the FIRST value for a duplicate key.
- `isCalibrated()` = `rt2ri.size() != 0` (0x6d2faa).
- `getRetentionTimeShift(rt)` (0x6d3566, vtbl[2]) and `getRetentionTimeCorrection()` both return `drt`. The shift is constant, not RT-dependent.

### 6.2 calculateRI(rt) (0x6d2b06, Observed) - piecewise linear, linear extrapolation

```
if rt2ri.size < 2: return -1.0
c = rt - drt
it = rt2ri.lower_bound(c)                  # first calibration RT >= c
if it->rt == c: return it->ri              # NB: dereferences end() if c > max RT (UB); control flow after that is unknown
if it == begin: ++it                       # below range -> use first segment
if it == end:   --it                       # above range -> use last segment
(rt1, ri1) = prev(it); (rt2, ri2) = it
return ri1 + (c - rt1) * (ri2 - ri1) / (rt2 - rt1)
```

This is the linear (van den Dool-Kratz) form only. No Kovats logarithmic form exists in the code: there is no other RICalculator implementation, and "ESTIMATED KOVATS RI" is only a UI label (Observed via symbols/strings). Extrapolation outside the calibration range uses the first or last segment.

### 6.3 calculateRT(ri) (0x6d2d34, Observed)

Inverse over `ri2rt` with the same segment selection: `rt = rt1 + (ri - ri1)*(rt2 - rt1)/(ri2 - ri1) + drt`. An exact key match returns `rt + drt`. With fewer than 2 pairs it returns -1.

### 6.4 Applying RI to compounds: calculateRI<int,float>(compounds) (0x3f98f4, Observed)

```
for comp in compounds:
    rt = comp.getRetentionTime()
    comp.setRetentionTimeCorrection(getRetentionTimeShift(rt))   # = drt
    comp.setRetentionIndex(calculateRI(rt))
```

Called at the end of every deconvolution (4.6). `getCorrectedRetentionTime()` (used in 5.6) is Inferred to be `RT - correction`.

### 6.5 RT shift correction: applyRetentionTimeCorrection<int,float>(comps, lib, ref, cutoff) (0x3f9e9e, Observed)

This is the "advanced RT shift correction". One reference library compound is located in the run, and the whole calibration is shifted so that the reference compound lands on its library RI.

```
win = 3 * LS_RI_DIFF; best = null; best_score = 0
for comp in comps:
    if comp.RI == -1: a = comp.RT/1000; b = ref.RT/1000  else: a = comp.RI; b = ref.RI
    if |a - b| < win:
        for h in lib.getLibraryHits(comp, (float) LS_RESOLUTION) while h.overall >= cutoff:
            if h.compound == ref and h.overall > best_score: best_score = h.overall; best = comp
if best == null: drt = 0
elif best.RI == -1: drt = best.RT - ref.RT
else: drt = calculateRT(best.RI) - calculateRT(ref.RI)      # RICalculator vtbl[1], Inferred = calculateRT
return drt
```

Callers then run `calculateRI(comps)` again:
- `RICalibrationWizard::RISetter::applyRICorrection` (0x3f5b32) uses its stored ref_comp and cutoff.
- `CalibrationSetter::applyRICorrection` (0x494f40) uses the reference compound from a combo box and the cutoff from a spin box.

Quirk (Inferred): `calculateRT` adds the OLD `drt`. A second application therefore computes the new drt relative to the old one, and does not accumulate it. The result is exact only when the old drt is 0.

### 6.6 Building the calibration: CalibrationPage::setDataToChromatogram (0x3eff6e, objdump; Observed)

```
ric = RICalculator("not available")
for row i in table:
    rt = toDouble(item(i, 1)); ri = toDouble(item(i, 2))
    if rt >= 0: ric.addCalibrationPair(MainWindow::mapXUnitToTime(rt), ri)   # display unit -> internal time
scan.setRICalculator(ric); scan.writeIndexFile()
```

`mapXUnitToTime` (0x3b5d6e) applies `TimeMapper::mapUnitToTime` only when a TIC plot exists.

### 6.7 Unknowns (RI)

- Wizard flow (strings, Observed): intro; chromatogram selection (file, "Recalibrate", "Use Retention time for chromatogram alignment (not recommended)" with a warning); "RI Calibration Mix" (library file, mix of library compounds, add/remove); "Calibration Table" ("Adjust the retention time to retention index mapping here"); then `RISetter` per chromatogram, with optional "advanced" reference shift (`CalibrationTablePage::performAdvancedCalibration/getLib/getCutoff`). `CalibrationSetter` (F9 dialog) applies the same to a list of chromatograms.

### 6.8 Mix assignment: CalibrationPage::matchCompounds (0x3eef58, objdump; Observed)

```
comps  = multimap{ comp.vslot2 (RT) -> hits }  for each detected component (library search without RI)
mix    = map{ c.vslot3 (RI) -> c }             for each mix compound (the set in the page)
found  = {}; it = begin(comps)
for (_, m) in mix (ascending RI):
    j = it; matched = false
    while j != end and not matched:
        (rt, hits) = *j++
        for h in hits: if h.overall >= cutoff (page field +0x68) and h.compound == m: found[m] = rt; matched = true
    if matched: it = j          # the next substance starts after this component
missing = mix - found
fillTable(found, missing)
```

`fillTable` (0x3ef57c) writes three columns (name, RT in display unit, RI) and sorts by RI. Found rows: RT from the component. Missing rows: RT = -1, red background, the RT cell edited through a `ComboDelegate` listing the RTs of all detected components. `setDataToChromatogram` (6.6) skips rows with RT < 0, so unmatched substances are left out of the calibration. The RT-only checkbox skips both pages: no `RICalculator` is set and the library search uses RT/1000 instead of RI (5.6).
- TimeMapper units, and the internal RT unit (ms is Inferred from the `/1000` fallbacks).
- RICalculator ctor default of `drt` (Inferred 0). The stream format (fromStream 0x6d315a) was not read.
- `recalibrateChromatograms` (0x3f52ec) was not read.

## 7. Quantification ion selection (Observed); other quantification NOT ANALYZED

### 7.1 Compound::autoQuantificationIonAssignment(n) (0x6109b0)

Calls slot 16 with `excl = QUANT_EXCL_IONS`, `min_dist = (float) QUANT_MIN_DIST_IONS` and `quality_cutoff = (float) QUANT_QUALITY_CUTOFF`.

### 7.2 Compound::autoQuantificationIonAssignment(excl, min_dist, n, quality_cutoff) (0x610442, Ghidra + objdump)

```
peaks = getPeakEssentials()                  # the PeakEssential copies stored in 4.6
if peaks empty: return {}
order = peaks sorted by integral DESC        # list::sort with peakEssentialCompare (0x5fa2c6: b.integral > a.integral), iterated in reverse
Q = set<float>(); count = 0
for p in order while count < n:
    if p.quality < quality_cutoff
       and not any(lo <= p.ion <= hi for (lo, hi) in excl)
       and all(|q - p.ion| >= min_dist for q in Q):
        Q.insert(p.ion); count += 1
if count < n:
    print "WARNING: Was not able to assign n quantification ions ... Using the first n with the highest intensity"
    for p in order while count < n: Q.insert(p.ion); count += 1     # restarts at the top, no filters
return Q
```

- Quality is "lower is better", consistent with the `quality <= 10` model-peak rule in 4.6.
- Fallback quirk (Observed): the fallback restarts from the highest-integral peak. Re-inserting an ion that is already selected still increments `count`, so the result can hold fewer than n ions.
- Duplicate ions from the duplicated peaks in 4.5 behave the same way.
- Ties in integral are ordered by list position (list::sort is stable, then reversed).

### 7.3 QuantIonChecker::detectQuantIons (0x5b70aa)

For every mapped compound `c`:

```
cc = c.getCentroidCompound()
added = cc.autoQuantificationIonAssignment(this.ion_number)
delete cc
if added is empty: drop the entry
else: comp2ions[c] = added
```

### 7.4 Unknowns (quantification)

- `LibraryCompound` overrides of the assignment (0x61db1c, 0x61e19a).
- `getCentroidCompound` (whether it keeps PeakEssentials).
- `QuantIonChecker::detectCompounds2Process`.
- The source of `ion_number` (Inferred QUANT_ION_NUMBER).
- The quantity reported per ion (Inferred: the PeakEssential integral from 4.9).
- MID: Section 7g (normalisation: Section 7c, statistics and PCA: Section 7f). Settings in play: QUANT_ION_NUMBER, QUANT_MIN_DIST_IONS, QUANT_QUALITY_CUTOFF, QUANT_EXCL_IONS. tools::Statistics::collectedMean (syms 0x6c3a2c), collectedVariance (0x6c383a / 0x6c3640) and polynomialRegression (0x6c3292) are entry points only.

## 7b. Natural isotope calculator (Observed)

`tools::NaturalMIDCalculator` (0x6d92bc) reads a table of `symbol, nominal mass, exact mass, percent` rows (`:/isotopes.csv`, embedded in the binary; extracted copy: `tools/isotopes/old-isotopes.tsv`). `calculateMID(formula, n)` (0x6d941e) returns `n` values: the product of one polynomial per atom, truncated, with index = nominal mass minus the sum of the lightest nominal masses (isotope gaps stay 0). `calculateExactMass` (0x6d9c6a) sums the lightest isotope masses. The formula grammar has no parentheses and no spaces. The table has typos (e.g. Sn 112 mass 111.05), so the rewrite uses IUPAC/NIST values instead (for about 30 common elements; the old table has 79, so rarer ones such as Ni, Cr, Sn, As, Pt still need rows); `crates/md-core/src/isotope.rs` is tested against the old output with the old values. Driver: `tools/isotopes/old_mid.cpp` (`LD_PRELOAD`, same idea as `tools/golden/`).

## 7c. Normalisation of the batch table (Observed; ported in `crates/md-core/src/normalize.rs`)

The values are the per-sample summed quant ion integrals `S[c][n]` of compound `c` in sample `n` (`CompoundSeries::getSummedIntegral`). The GUI offers four modes; z-score is a flag on top of any of them (`BatchQuantificationTreeModel::getTableData(zscore, exclude_norm_row)`, 0x56527a). The data() cells (0x565a42) use the same divisors for the per-ion rows.

```
none:      v = S[c][n]
internal:  v = S[c][n] / S[std][n]            # std = BatchQuantificationTreeModel::norm_series (normalizeData, 0x5657f8); its own row becomes 1; no zero guard
sum:       v = S[c][n] / T[n]                 # T[n] = ResultSeries::getTotalSampleIntensity(n) (0x734208) = sum over ALL compounds of getSummedIntegral(n), the standard included; normalizeDataBySum 0x5658a8 clears norm_series (the two modes exclude each other)
pool:      v = S[c][n]                        if n is a pool sample
           v = S[c][n] / ref                  else   (PoolNormalizer::getNormalizedSummedIntegral 0x72e194, per-ion variant 0x72de96)
             ref = mean | median of S[c][p] over p in pool_mapping[n]
z-score:   z = (v - mean(row)) / sqrt(variance(row))    # row = one compound over all samples, after the mode above
```

- **Pool samples:** the replicate group named `$POOL` (`BatchQuantificationResultWidget::normalizationByPool` 0x554a14, `getGroupIndices("$POOL")`). The widget also reads the acquisition date of every sample file into `PoolNormalizer::times` (long long); the value type is a QDateTime read, the unit is not verified. The pool samples keep their raw values.
- **Pool mapping** (`updatePoolMapping` 0x72e52c): for every sample `n`, all pool samples go into a `multimap<|times[n] - times[p]|, p>`; the first `count` entries in key order (equal distance: the lower pool index first, multimap insertion order) form `pool_mapping[n]`. `count` is a PoolNormalizer field; the dialog default is not read (the rewrite uses 2). A pool sample is in the mapping too (its own nearest pool is itself).
- **Pool mapping export** (`BatchQuantificationResultWidget::exportPoolMapping` 0x5552a0): a text file, `.txt` forced on the name. Content: `Mapping`, newline, `Sample<TAB>Pools`, two newlines, then one line per sample in the table's order: the sample name, and for each pool in `pool_mapping[n]` a TAB and the pool's name, newline. The rewrite writes the runs the table shows (`normalize::pool_mapping_text`).
- **Mean/median** (`mode` MEAN or MEDIAN): `gsl_stats_mean`, and `gsl_sort` + `gsl_stats_median_from_sorted_data` (the middle value of an odd count, else the mean of the two middle values). Missing values are stored as 0 in the old table and count in these statistics.
- **Z-score:** `tools::Filter::zScoreTransformation<double>` (0x55a1aa, element formula confirmed in the disassembly: `(x - mean) / dev`), with `Statistics::variance` (0x7c2efe) = `sum((x - mean)^2) / (n - 1)` (sample variance, N-1).
- **Division by zero:** nothing is guarded, results are inf/NaN. The rewrite returns `None` instead (the GUI shows a dash).

Deviations of the rewrite (the table has `None` for "not found" instead of 0):

- A run without a value stays `None` after every mode.
- For `total`, a missing value adds nothing (same as the old 0). Substances dropped by the reproducibility filter are not in the table, so they are not in the sum either (not checked whether the old program sums before its filter).
- Missing values of a pool run are skipped in the mean/median, and missing values are skipped in the z-score statistics (the old 0 would drag the mean down). `z_score` needs two values and a non-zero spread, else the row is `None`.
- Acquisition time: the app uses the file date; the old program reads the date stored in the file.
- Not ported: the old program scales each ion row of the tree view; the table here has one summed row per substance.

### 7c.1 Unknowns (normalisation)

- Dialog defaults for the pool `count` and MEAN/MEDIAN mode, and the unit of `times`.
- `ResultSeries::getTotalSampleIntensity` sums `getSummedIntegral`, which would apply the PoolNormalizer when one is set; the GUI never combines pool and sum.

## 7d. Spectrum display normalisation (Observed)

`SpectrumPlot::actNormalizationChanged` (nm 0x4156d2; Ghidra addresses are nm + 0x100000) maps four toolbar actions to `Normalization` = `NONE`, `SUM`, `SSUM`, `BP` and calls `setNormalizationMode`, which re-plots `getNormalizedData` (0x4177de) for the profile and the centroid spectrum. `getNormalizedData` divides every integer intensity by one `norm` (double division, no rounding):

| Mode | `norm` |
|------|--------|
| `NONE` | 1 |
| `BP` | `getBasePeakIntensity()` (0x601d22), the largest intensity |
| `SSUM` | `getSpectrumVectorLength()` (0x601c32) = `sqrt(Σ i²)` (decompiled) |
| `SUM` | `getTotalSignal()` (0x601b6e) = Σ i (name only, not decompiled) |

The rewrite calls `SSUM` "Unit" (unit vector) and `SUM` "Relative" (fractions of the total signal), `md_core::spectrum::normalize`. A zero norm is not guarded in the old code (division by zero); the rewrite leaves such a spectrum unchanged. The toolbar label strings were not read, so the names Unit/Relative are inferred from the issue. `SplitSpectrumPlot` and `LibrarySpectrumPlot` reuse the same mode.

## 7f. Result table of a run and Find Compound (Observed; ported in `crates/md-core/src/results.rs`)

Observed with objdump on `DataTableModel::getData` (0x4c88a2), `CompoundExporter::run` (0x4b5bb8) and `CompoundMapper::findCompound` (0x4b6ca0).

- Columns: Compound Name, RT, RI, Peak Area, Overall Score, Spec. Score, RI Score, Quant. Ions. A component without hit (or a hit below the model's cutoff) shows `No match: ` plus its own name, the scores as -1 and the area as -1.0. RT goes through `TimeMapper::mapTimeToUnit`.
- Peak Area (`getIntegral`): sum over the hit library compound's quant ions of the component's integral of that ion (tolerance float constant at 0xa824e8, value not read). If the library compound has no quant ions, the component is asked for three ions (slot 0x88, argument 3; the meaning is Inferred, ported as the three largest peaks).
- Quant. Ions: the library compound's quant ions as `QTextStream` floats, fixed notation, precision 4, joined by ` / `. The old code truncates two characters from the end, which leaves a trailing space. Empty when the set is empty.
- `CompoundExporter` writes the model generically: the header row, then every cell's display text. The separator is a TAB (0xa810ab), the line end `\n`, no quoting, despite the `*.csv` file filter. The new export uses comma and RFC 4180 quoting and leaves missing values empty instead of -1.
- Find Compound: for each component's list of hits, keep a hit whose library compound has the same name (string compare) as the chosen one and an overall score of at least `LS_CUTOFF`; list (score, component), sorted descending (`qSort` with `qGreater`).

## 7g. Component panel and split spectrum plot (ported in `crates/md-core/src/detail.rs`)

Observed with objdump on the component panel labels, `Compound::getAveragePeakQuality`, `getAverageSignalToNoise`, `getPeakEssentials` and `SplitSpectrumPlot::setSplitData` (0x41d90a).

- Peak dis.indx (peak discrepancy index): `getAveragePeakQuality`, the mean of `PeakEssential::getQuality` over the compound's `map<float mass, PeakEssential*>`. 0 is best. An empty map gives -1.0. The same value is needed for the quality filter of #58.
- S/N: `getAverageSignalToNoise`, the mean of the signal-to-noise over the same map; 0 for an empty map.
- Peak sharp.: the mean of `getSharpness` over `getPeakEssentials()`. The tooltip says "highest", but the code takes the mean.
- The map holds one peak per mass. Two peaks of the same mass keep the first, which is the sharpest (§4.5). The rewrite does the same (`component_detail`).
- The old panel lists at most 10 hits; the rewrite keeps the 15 of the hit table.
- Split plot: `setSplitData` builds two curves. The first one is the found spectrum as is, the second one is the library spectrum with every intensity negated (`xorpd` with the sign mask), so the library hangs below the axis. Both are set with `QwtPlotCurve::setSamples`. Observed.
- Not decoded (shortcuts in the rewrite): how the old difference mode computes its data, the `SpectrumMode` values, and the layout of `exportCSV` / `exportToClipBoard`. The rewrite subtracts the library from the found spectrum after the same normalisation (without one the base peak, because the found heights are counts and the library is scaled to 0..999; found minus library, on the union of both nominal masses) and exports both spectra and their difference as CSV; the clipboard gets the same text.
- Ion counts: the rewrite reads centroid spectra only, so "ions" and "centroid ions" are the same number, and there is no profile switch (architecture.md).

## 7e. Adding a found compound to a library, copying RTs (Observed; ported in `crates/md-core/src/library_add.rs`)

Read from the disassembly of `DialogAddLib::on_okButton_clicked` (0x4ca1be), `Ui_DialogAddLib::retranslateUi`, `LibrarySearch::addCompoundSpectrum` (0x627f88) and `MainWindow::copyRT` (0x3cacc2).

**Add dialog.** Input is the selected component (`Compound`), the target `LibrarySearch`, a name (editable combo box of the library names) and a comment. Without an open library: error "No compound library opened. Changes can not be saved." The name must not be empty. The dialog copies the compound, sets the name to the combo text and, when the comment is not empty, adds the feature `COMMENT`. If the library does not contain the name, the compound is added (`addCompound(c, "no source")`). Otherwise one of four radio buttons decides, each after a yes/no question:

| Button | Text | Action |
|---|---|---|
| `rb_replace` | Replace entry | Removes the entry and adds the found compound (question "replace the existing reference compound?"). |
| `rb_modify` | Modify entry | `addCompoundSpectrum(c, name, source)`: merges the spectrum into the entry (combined spectrum); not ported. |
| `rb_keep` | Replace RT/Spec | A new `LibraryCompound(c, ...)` takes the auto quant ions, the quant ions, the features and the formulas of the old entry and replaces it: spectrum, RT, RI new; ions and comments old. |
| `rb_quant_ions` | Exchange quant. ions | Sets only the auto-detected quant ions of the entry to those of the found compound. |

After a change the library is written with `toDisk()` ("Error saving compound library" on failure). The auto quant ions of the found compound come from `autoQuantificationIonAssignment(n)` (§7.1).

**Copy retention times** (File menu, "Copy Retention Times"): asks for a target library (`.lbr`), maps the components of the run to it with `CompoundMapper::mapCompsToLib` (same search and one-to-one assignment as §5), and shows a table with the columns Compound, RT, Lib RT, RI, Lib RI, Score. Entries with several candidates are set in italics. The write-back after confirmation is in a slot that was not read; the rewrite replaces RT and RI of the assigned entries.

Deviations and gaps of the rewrite:

- "Modify entry" (merge of spectra) is not ported; the combined-spectrum arithmetic (`LibraryCompound::addSpectrum`, virtual call from 0x61cde2) was not read.
- With "Replace RT/Spec" the old code copies the old features over the new ones, which drops a comment typed in the dialog; the rewrite appends the typed comment (new line) to the old one.
- The old "Replace RT/Spec" also copies the formulas of the old entry; `LibraryCompound` has no formula field, so they are not copied. "Replace entry" removes the entry and adds the new one at the end; the rewrite replaces it in place. Like the old dialog, the entry gets the source "no source"; the run name is only in the suggested entry name.
- A found compound without RI (run not calibrated) is stored with RI -1 (no RI); copying RTs then leaves the RI of the entry alone.
- The result is an entry in the rewrite's `.mdlib` (saved through the library editor), not an `.lbr`.
- Source string "no source" and the `LibraryCompound(Compound, source)` constructor defaults for the combined spectrum are Inferred from the MSL importer.

## 7f. Biomarker statistics, box plots and PCA (Observed; ported in `crates/md-core/src/stats.rs`)

Addresses are `nm` addresses of the ELF (Ghidra = nm + 0x100000). Found by disassembling `BiomarkerWidget::exportData` (0x4882e6; `createGUI` 0x48a044 repeats the same loop for the table), `BoxPlot::setData` (0x48d42c), `BatchQuantificationResultWidget::showPCA` (0x4557e0) and `tools::PCA` (0x6d7bc2).

**Two-group table** (`BiomarkerDialog` picks two lists of replicate groups; `BiomarkerWidget(series, group1, group2, independent)`; the flag is the byte at `this+0x60`). For every substance:

```
independent: x = getSummedGroupIntegrals(g, compound) for every g in group 1   # one value per run
             y = the same for group 2
else:        x = getGroupMean(g, compound, true) for every g in group 1        # one value per replicate group
             y = ... group 2
hits1, hits2 = sum over the groups of getGroupHits(g, compound)                # runs of the group with a hit
mean1 = alglib::samplemean(x), mean2 = alglib::samplemean(y), ratio = mean1 / mean2
p = alglib::unequalvariancettest(x, nx, y, ny, &both, &left, &right) -> both   # two-sided
```

- The only test is **Welch's t-test** (unequal variance, two-sided); `studentttest1/2` are linked but never called from the dialog. alglib 3 (0x744ef0): `nx < 2` or `ny < 2` gives p = 1; a zero sample variance in either group gives p = 1; else `t = (mx - my) / sqrt(vx/nx + vy/ny)`, `df = (vx/nx + vy/ny)^2 / ((vx/nx)^2/(nx-1) + (vy/ny)^2/(ny-1))` (not rounded) and `p = I_{df/(df+t^2)}(df/2, 1/2)` with the regularised incomplete beta function (`incompletebeta` 0x7b7800), which is the same p-value as R's `t.test(x, y)`.
- There is **no paired test and no multiple-testing correction** in the old program. "Independent samples = no" means a Welch test over the replicate group means. The rewrite adds a paired t-test (group k of side 1 against group k of side 2 in the order given, run by run; same number of groups and of runs per pair, else no p-value; `t.test(paired = TRUE)`) because issue #24 asks for it.
- Missing values are 0 in the old table, so they are part of the means and of the test. The ratio is not guarded (division by zero gives inf/NaN).
- CSV (`exportData`): tab separated, header `Metabolite  Identification Score  Hits 1  Hits 2  RT  RI  Mean G1  Mean G2  Ratio  P-value`, then the values of group 1 and of group 2. "Identification Score" is the reference compound's `lib_score` feature. Further sub-features: paste a sample list to fill a group (`autoSelect`, exact match of each item with a group name), filter by p-value, choose groups and their order for the box plot.

**Box plot** (`on_b_boxplot_clicked` 0x421842, `BoxPlot::setData`): per selected substance and per selected group the `getSummedGroupIntegrals` values (all runs of the group, sorted). `QCPStatisticalBox::setData(key, min, q1, median, q3, max)` with `min`/`max` the smallest and largest value (no outlier rule), `median = gsl_stats_median_from_sorted_data`, `q1`/`q3 = gsl_stats_quantile_from_sorted_data(0.25 / 0.75)` (linear interpolation, `h = (n-1)p`, R `quantile(type = 7)`). Boxes sit at x = 0, 1, 2, ... in the chosen order.

**PCA** (`showPCA` and `tools::PCA`):

```
rows = getTableData(false, true) of the shown table    # substance x run (as normalised on screen)
runs = all runs except the replicate group "$POOL"
a row is dropped if any value is inf or NaN; else it is z-scored over the runs
        (tools::Filter::zScoreTransformation, N-1 variance; a constant row gives NaN and is dropped)
PCA(dims = the kept rows, samples = the runs):
    centre every dimension (gsl_stats_mean subtracted; a no-op after z-scoring)
    C = covariance matrix of the dimensions (gsl_stats_covariance, N-1)         # = correlation matrix here
    gsl_eigen_symmv(C) ; eigenvectors * -1 ; gsl_eigen_symmv_sort(GSL_EIGEN_SORT_VAL_DESC)
getTransformedData(k)   = centred data * first k eigenvectors                    # scores, runs x k
getPrincipalComponent(i)= eigenvector i                                          # loadings, one value per substance
getProportionOfVariance = eigenvalue_i / sum(eigenvalues)
```

It throws `PCAException` for an empty input or dimensions of unequal length. `PCAPlot` has combo boxes for the x and y component (`comboPCAChanged`), a score plot, a "dimension impact" (loadings) plot and a histogram plus cumulative histogram of the proportion of variance.

Deviations of the rewrite:

- The eigenvectors are computed through an SVD of the centred runs x substances matrix (`nalgebra`). It has the same components as the covariance eigenproblem (`prcomp(t(Z), scale = FALSE)` in R) but costs `O(runs^2 x substances)` instead of `O(substances^3)`. Components with a numerically zero eigenvalue (beyond runs - 1) are kept with proportion 0.
- The old program flips every eigenvector (`* -1`), GSL's own sign being arbitrary. Signs are arbitrary in any implementation; the rewrite fixes the sign so that the largest absolute loading of each component is positive. The tests compare absolute values or the sign-invariant product.
- A missing value is `None` in the rewrite and 0 in the old table; PCA and the two-group test treat it as 0 by default (`missing_as_zero`), or skip it when that option is off. On a z-scored table a missing value is always skipped (0 would mean "average" there). Values that the normalisation could not compute are also `None`.
- Not ported: PDF/clipboard export of the plots (PNG only).

### 7f.1 Unknowns (statistics)

- Whether `BiomarkerWidget::createGUI` differs from `exportData` besides the table widget (the same alglib calls appear at 0x48a94a-0x48ab20).
- The default p-value filter and the dialog defaults for the independent flag.

## 7g. MID wizard (Observed; ported in `crates/md-core/src/mid.rs` and `mid/`)

Addresses are nm addresses (Ghidra = nm + 0x100000). Feature F16.

### 7g.1 Inputs and method
`MethodPage::getSelectedMethod` (0x4f4b00): 0 uncorrected, 1 corrected by unlabelled reference runs (default), 2 by library spectrum, 3 by sum formula. `getM1Correction` (0x4f4b6c) is the abundance spin box when the M1 check box is on, else 0; the wizard divides it by 100. `getTracer` (0x4f4bbc) is the element of the tracer combo box. `MIDWizard::accept` (0x4f8354) makes the groups `GL` (labelled runs) and, for method 1 only, `GUL` (unlabelled runs). `createMIDLibrary` (0x4fa038) keeps the library substances with fragments and makes every fragment ion `lo..=hi` a quant ion, so the batch rescan (`ResultSeries::preciseSeriesCompletion`, Section 7) fills ions that were not deconvoluted. The settings page is the batch settings widget in targeted mode. Port: `mid::mid_library`, then the normal batch scan with rescan.

A fragment is a key `(a, b)` of the library's sum formulas (`LibraryCompound::getFragments` 0x61f428, map order); `lo = min`, `hi = max`, ions `lo..=hi`. `getFormula` and `getMIDSize` (default 0 = one share per ion) belong to the same key. Port: `mid::fragments`.

### 7g.2 Group methods 0 and 1: `Calculate::operator()` (0x4fa57c)
Per substance and fragment, one row of integrals (`CompoundSeries::getIntegral` 0x6319c8, 0 for a missing ion) per run of `GL` (and `GUL`). Rows whose integrals are all 0 are not used. Method 0: `MIDSolver(rows)` (0x6e00cc): mean of the normalised rows, CI = 1.96 · sd / √n (sample sd), CI -1 for one row, R² -1. Method 1: `MIDSolver(GL rows, GUL rows, m1)` (0x6dfe5a): every pair (reference, labelled) adds one equation per ion; column k is the reference shifted by k, with the M+0/M+1 entries corrected by `cf = k · m1` as `shift = r0 · cf / ((r1/r0 + 1) − cf)` (diagonal `r0 + shift`, sub-diagonal `r1 − shift`, both clamped at 0). Errors (`GCMSScanException`) become the feature `ERROR`, "Error calculating reference spectrum for <b>name</b> (lo-hi): message"; `midCalcFinished` (0x4f9534) shows them in one WARNING box. Port: `mid::determine` with `Method::Uncorrected` / `Method::Reference`, solver `correct_by_reference`, `replicate_mean`.

`Calculate` (ctor 0x4fbb28) holds the `ReplicateResultSeries`, the method, M1, tracer and calculator; it reads the runs of `GL`/`GUL` with `getGroupIndices` (sample names of the group map at +0x150 → indices into the sample names at +0x80) and the integrals with `getCompoundSeries` (map at +0x20). It returns a `MIDCompound` (0x4f6366, a copy of the library substance with one `MIDSolver` per fragment key at +0x290). Only `MIDCalculatorException` is caught. With no usable reference row or no usable labelled row, method 1 gets the solver's empty result (MID [-1], CI [-1], R² -1, no error). With no usable labelled row, method 0 throws `std::out_of_range` (`rows.at(0)` in 0x6e00cc), which leaves `Calculate`; the port gives no MID there, and its export writes the sentinel row of the empty result (§7g.7) where the old program writes none.

The old program never shows these results: `midCalcFinished` builds the `ReplicateMIDData` only for methods 2 and 3, and `MainWindow::midWizClosed` (0x3ca9a8) opens the result window (`MIDSingleResultWidget`) only when it exists; for methods 0 and 1 it shows "No compounds found during MID determination. Check library and RI window." `MIDResultWidget` (0x4e9104), the group result window, has no caller. The port shows the group results (deviation).

### 7g.3 Per-run methods 2 and 3: `MIDCompoundSeries` (ctor 0x6edb74, `determineMIDs` 0x6ee5e0)
Every run of the batch gets its own MID (`Calculate_Individual` 0x4fb7fa; flag = method 2). Per fragment and run:
- Method 2: reference row = `getNormalizedIntensity(lo + k)` of the library spectrum (normalised to 1000, exact float key, 0 when absent), correction = M1/M0 of the tracer's natural pattern (`calculateMID(tracer, 2)`), `MIDSolver(row, reference, correction)` (0x6dfd1e).
- Method 3: `MIDSolver(row, formula, calc, tracer, size)` (0x6dff04 → 0x6e1896): column k is the normalised natural pattern of the formula with k tracer atoms removed (none once the count is 0), fit by least squares.

Results are keyed by ion `lo + n` (MID, CI, R² per run). An empty tracer or formula throws "Invalid formula string"; the whole substance is dropped and the wizard reports "Error determining MIDs.". A run without the substance has `getIntegral` 0 and gives NaN in the old program.

### 7g.4 Least squares (`MIDSolver::solve` 0x6e1148 / 0x6e0656)
`gsl_multifit_linear`; R² = 1 − χ²/TSS; when rows ≠ columns, CI_k = √(COV_kk · χ²/(n−p)) · t₀.₉(n−p), else -1. Empty input gives the sentinels MID [-1], CI [-1], R² -1. Port: `mid::fit` (Householder QR, t quantile by bisection).

### 7g.5 "Export MIDs" of the batch table (`BatchQuantificationResultWidget::exportMIDs` 0x452132)
Uses the feature `FRAG`, parsed by `parseFragmentString` (0x45887e): regex `\((\d+)/(\d+)\)`, (first ion, number of ions) into a `QMultiMap` (sorted by ion; a repeated ion lists its latest entry first). The file has three tables, sample names = file names without `.bin`:
1. "Retention Times (Min)": one row per substance with `FRAG`, RT in minutes or -1.
2. "Total ion counts": rows `Name_ion`, the sum of the integrals of `ion..ion+size-1` per run. A substance without `FRAG` gives "Warning: No fragment (FRAG) defined for compound X"; a `FRAG` that does not parse gives `Error: Can not parse "…". SYNTAX(e.g): (174/8)(261/5)`.
3. "Mass Isotopomer Abundances": rows `Name_m (Mk)`, integral / fragment total, a blank line after each fragment.
Port: `mid::frag_export`, `mid::parse_fragments`.

### 7g.6 Verification and deviations
- Solver: 31 synthetic cases of `tools::MIDSolver` called through `LD_PRELOAD`; test `matches_old_program` in `mid.rs` (max error 4e-15).
- Methods 2 and 3 end to end: `tools/isotopes/old_midseries.cpp` drives `MIDCompoundSeries` on synthetic integrals (`make-midseries.sh`, `gen_midseries.py`); `crates/md-core/tests/mid_series.rs` matches all 90 values (MID, CI, R²; several fragments, reversed fragment keys, MID size, tracers C and N, no tracer, runs without integrals), with the library spectrum read from the old `lib_aa.lbr`.
- Methods 0 and 1 end to end: `tools/isotopes/old_midgroups.cpp` builds a `ReplicateResultSeries` with the groups `GL`/`GUL`, runs `Calculate` and the export of the (otherwise unused) `MIDResultWidget` (`make-midseries.sh`, cases from `gen_midseries.py groups`); `tests/mid_series.rs` matches all 45 values (MID, CI, R²; several fragments, reversed keys, missing ions, runs without integrals, one labelled row, square method 1 fit, M1 correction, no reference row, no labelled row) and the 11 exports byte for byte.
- Deviations: a run without the substance or with only zero integrals gives no MID (`None`) instead of NaN (the exports print `nan` as the old program); the FRAG export keys fragment totals per row instead of per substance id; RT minutes are `rt_s / 60`. Fragments wider than 1000 ions (`MAX_FRAGMENT_IONS`, a library key span or a `FRAG` size) are left out; the old program loops over all of them. The app's MID window looks up the RT shift reference in the fragment-only scan library and scans without shift when the reference has no fragment.

### 7g.7 Result window: group statistics and text exports
Numbers are written by `QTextStream << double`: `%g` with 6 significant digits, `nan`, `inf` (port: `mid::qt_g`, also used by the FRAG export of 7g.5). Run names are `QFileInfo::baseName` (file name up to the first dot).

- Statistics (`ReplicateMIDData`, built once in `MIDWizard::midCalcFinished` 0x4f9534 for the per-run methods; one group `GL` with all runs). `getGroupMIAbundances` (0x6f6052 mean, 0x6f64ee standard error) collects `MIDCompoundSeries+0x90[i].at(run)` (one vector per isotopomer, fragments in order); every caller passes the filter flag, which drops values that are 0, NaN or infinite. Mean = `gsl_stats_mean`, standard error = `gsl_stats_sd / √n`. No value gives mean 0, fewer than two give a NaN standard error. Port: `mid::group_stat`.
- `MIDSingleResultWidget::outputMIDDataToStream` (0x4ef548, per-run methods): header `\t<run>` per run, `\t`, `\tCI (<run>)` per run, `\t`, `\tR2 (<run>)` per run. Per substance `Name;RI;RT` (library spectrum, RT in ms), per isotopomer a line break before `M0`, then `ion (Mk)`, the MIDs, `\t`, the CIs, `\t`, the R² of each run; two line breaks after the substance. Isotopomers per fragment: one per ion, at most the MID size when it is set. A run without integrals is `nan` for MID and R²; its CI is `nan` for the sum formula and -1 for the library spectrum (method 2 has no CIs). Port: `Mids::run_text`.
- `outputMIDStatisticsDataToStream` (0x4f0090): header `\tGL\t\tStdErr (GL)\tMetabolite name\tFragment ion\tMass isotopomer`, the same substance lines, rows `ion (Mk)\tmean\t\tstderr\tName\tlo\tMk`. Port: `Mids::stats_text`.
- `MIDResultWidget::actExportTriggered` (0x4e9212, group methods): no header; per substance and fragment key the rows `Name_ion (Mk)\tMID\tCI\tR2` (CI and R² -1 where there are none) and an empty line after each fragment. The rows stop at the solver's MID count, so the empty result is one row `M0` with -1, -1, -1. Port: `Mids::group_text` (a fragment without MID writes that row). Deviation: the old loop runs from the key's first ion to its last, so a reversed key `(b, a)` writes no rows; the port writes them.
- Verified: `tools/isotopes/old_midstats.cpp` (run by `make-midseries.sh`) builds the `ReplicateMIDData` and the result widget on the same cases and prints both exports; `tests/mid_series.rs` matches all 34 statistics and both exports byte for byte; the group methods and their export are verified in 7g.6.

## 8. Open gaps (ranked by impact on reproducing results)

1. Default setting values. The storage format is known (Section 1); the values are not in code read so far. Find `fromDisk` callers and any shipped settings file, or the GUI defaults. This blocks golden-output comparison for every section.
2. (resolved) `tools::Filter::matched2ndDerivativeGaussianFilter` edge behaviour: the ends hold the nearest inner value (Section 3.3).
3. Spectrum binning and centroiding (`getBinnedSpectrum`, `getCentroidCompound`). Library scores depend on them when bin widths differ.
4. DefaultPeakDetector details:
   - the C=0.2 tail rule (value read from rodata, not verified against output)
   - the `i` update after a rejected peak (possible off-by-one)
   - the identity of vtbl[4]/[5]
5. MatchedPeakDetector: resolved (Section 3.3).
6. (done, 6.8)
7. Pooled statistics formula (collectedMean/collectedVariance; N-1 denominator Observed) and the meaning of the `n` field in BaselineStatistics.
8. The source of the noise factor that GCMSAnalyzer passes to setNoiseFactor (Section 2.3).
9. (done, 7g.6) MID group methods 0 and 1 and their export against the old binary. Normalisation (7c), statistics and PCA (7f) are done.

## 9. Ranked next steps

- Locate `GCMSSettings::fromDisk` callers and any default settings file in the install.
- Decompile `Spectrum<int,float>::getBinnedSpectrum` / `getCentroidCompound`.
- Produce golden output from the old binary on one sample (peaks, components, spectra, hits, RI) and diff the Rust port against it. Reproduce the quirks in 4.5, 5.3 and 7.2 first, then decide whether to fix them.

Status: Deconvolution complete (Section 4). Library scoring complete (Section 5; binning/centroid internals open). RI calibration complete (Section 6, mix assignment 6.8). Quantification ion selection complete (Section 7.1-7.3). Normalisation reconstructed (Section 7c). Statistics and PCA reconstructed (Section 7f). MID reconstructed (Section 7g).
