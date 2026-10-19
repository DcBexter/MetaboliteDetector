# MetaboliteDetector 3.3 (Linux x86-64) - algorithm notes

Status: partial. Noise/baseline, peak detection and peak metrics (ported in `src/md-core/src/peaks.rs`, deviations in `docs/architecture.md#peak-detection`), deconvolution (Section 4), library scoring (Section 5, ported in `src/md-core/src/search.rs`), RI calibration (Section 6) and quantification ion selection (Section 7) are reconstructed. MID, normalizers, statistics and PCA are NOT analyzed yet (see Section 8).

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
- Keys also seen in code but not in toStream: AN_MIN_HEIGHT_BP, AN_MODELED_VALUES_FOR_INTEGRATION, SQL_MINIMAL_NUMBER_EQUAL_FRAGMENTS, SQL_MAX_NUMBER_FRAGMENTS, LS_CUTOFF, LS_MASS_FILTER, LS_IDENTIFICATION_ALGORITHM. Their persistence is Unknown.
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

```
W = peak_width (ctor arg); thr = ctor threshold (int)
f = tools::Filter::matched2ndDerivativeGaussianFilter<int>(values, density)   # Observed call; density = 4*W is NOT confirmed
regions: begin at i where f[i-1] < thr < f[i]; end at first j > i with f[j] <= thr
accept region if (end - begin) > 2
createPeak(ion, begin, end, f[begin..end])
```

- Filter (syms 0x64a5fd int, 0x64966d double; Observed): for x in [-d, d], u = 4x/d, kernel k(x) = (1.73/sqrt(d)) * (1 - u^2) * exp(-u^2/2), a negative 2nd-derivative Gaussian with sigma = d/4. Taps outside the signal count as 0; the int variant truncates the result (cvttsd2si). Edge behaviour beyond that: Unknown.
- Ctor (0x5ee926) takes (size_t peak_width, int threshold) without defaults. Its only caller is CalibrationSetter::threadFinished (0x4947f4), with both values from spin boxes of the RI calibration dialog (Observed; pairing Inferred).
- The region rules above are Inferred.
- Whether createPeak gets the raw or the filtered vector is ambiguous. The decompile passes the filtered vector (Observed), so the peak area is taken from the filtered data (Inferred interpretation).

### 3.4 MatchedPeakDetector<int,float>::createPeak (MCP 0x6eed8e; syms 0x5eed8e)

- apex = argmax over [begin, end] (strict >) on the supplied values.
- Peak(ion, begin, end, apex, values, nf = 0.5). Observed.

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
- Batch quantification: `ResultSeries::scanIons` (0x638589..0x6385d7) builds an analyzer with a fresh `DefaultPeakDetector`, calls `setMinimumPeakNumberRequired(1)`, then `deconvoluteChromatogram(peaks, scan_count, start_scan, 2 * AN_DECONVOLUTION_WIDTH)` on a scan window.
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

Observed (tail at 0x5e8e99..0x5e8f15). The meaning of data type 2 is Unknown. Deconvolution passes 1.0 as the Compound "resolution" for type 2, else 0.0 (Section 4.6), so type 2 is most likely nominal-mass (unit resolution) data (Inferred).

### 4.2 reorderPeaks(scan_dist, ion_dist) (0x5eac2c) - align apexes of split masses

Used for non-type-2 (high-resolution) data, where one nominal fragment can be split over neighbouring m/z traces. With ion_dist = 0.05 it never merges unit-mass traces. Observed (Ghidra decompile, structure confirmed):

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

- `insertIonProfile(bin, bins, dist, p)` (0x5eb972, objdump): `prof = new IonProfile(p)`. Then `lo = (int)(bin - dist) < 0 ? 0 : bin - dist` and `hi = min(bin + dist, bins.size - 1)`. It pushes `prof` into `bins[k]` for k in [lo, hi). The upper bound is exclusive, so bin `hi` is not filled (Observed: `jae` at 0x5eba59).
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

- `tools::Filter::matched2ndDerivativeGaussianFilter<double>` edge handling (kernel: Section 3.3). This is needed to reproduce the threshold T = 5.0 exactly.
- Meaning of `getDataType() == 2`, and how `Compound(rt, ints, ions, res)` uses `res`.
- `scan->vtbl[3]` (scan count) and `vtbl[4]/[5]` (ion count / EIC) identities are Inferred from usage.
- Whether the duplicated bin peaks from 4.5 are intended. A port should reproduce them first to match golden output.

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
- `LS_IDENTIFICATION_ALGORITHM` and `LS_USE_SCALED_LIB` semantics.
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

- How the wizard proposes the alkane RT/RI table rows (CompoundDetector + library matching in RICalibrationWizard, startRICalibration 0x3e93b8). Not read.
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
- MID, normalizers, statistics tests and PCA are not analyzed. Settings in play: QUANT_ION_NUMBER, QUANT_MIN_DIST_IONS, QUANT_QUALITY_CUTOFF, QUANT_EXCL_IONS. tools::Statistics::collectedMean (syms 0x6c3a2c), collectedVariance (0x6c383a / 0x6c3640) and polynomialRegression (0x6c3292) are entry points only.

## 7b. Natural isotope calculator (Observed)

`tools::NaturalMIDCalculator` (0x6d92bc) reads a table of `symbol, nominal mass, exact mass, percent` rows (`:/isotopes.csv`, embedded in the binary; extracted copy: `tools/isotopes/old-isotopes.tsv`). `calculateMID(formula, n)` (0x6d941e) returns `n` values: the product of one polynomial per atom, truncated, with index = nominal mass minus the sum of the lightest nominal masses (isotope gaps stay 0). `calculateExactMass` (0x6d9c6a) sums the lightest isotope masses. The formula grammar has no parentheses and no spaces. The table has typos (e.g. Sn 112 mass 111.05), so the rewrite uses IUPAC/NIST values instead (for about 30 common elements; the old table has 79, so rarer ones such as Ni, Cr, Sn, As, Pt still need rows); `src/md-core/src/isotope.rs` is tested against the old output with the old values. Driver: `tools/isotopes/old_mid.cpp` (`LD_PRELOAD`, same idea as `tools/golden/`).

## 8. Open gaps (ranked by impact on reproducing results)

1. Default setting values. The storage format is known (Section 1); the values are not in code read so far. Find `fromDisk` callers and any shipped settings file, or the GUI defaults. This blocks golden-output comparison for every section.
2. `tools::Filter::matched2ndDerivativeGaussianFilter` edge behaviour (kernel known, Section 3.3). It sets the component boundaries in 4.5 together with the fixed threshold 5.0.
3. Spectrum binning and centroiding (`getBinnedSpectrum`, `getCentroidCompound`). Library scores depend on them when bin widths differ.
4. DefaultPeakDetector details:
   - the C=0.2 tail rule (value read from rodata, not verified against output)
   - the `i` update after a rejected peak (possible off-by-one)
   - the identity of vtbl[4]/[5]
5. MatchedPeakDetector: the density it passes to the filter, its region rules, and raw or filtered values in createPeak.
6. The RI wizard's automatic alkane assignment (6.7).
7. Pooled statistics formula (collectedMean/collectedVariance; N-1 denominator Observed) and the meaning of the `n` field in BaselineStatistics.
8. The source of the noise factor that GCMSAnalyzer passes to setNoiseFactor (Section 2.3).
9. MID, normalizers, statistics tests, PCA: not started.

## 9. Ranked next steps

- Locate `GCMSSettings::fromDisk` callers and any default settings file in the install.
- Decompile `Spectrum<int,float>::getBinnedSpectrum` / `getCentroidCompound`.
- Produce golden output from the old binary on one sample (peaks, components, spectra, hits, RI) and diff the Rust port against it. Reproduce the quirks in 4.5, 5.3 and 7.2 first, then decide whether to fix them.

Status: Deconvolution complete (Section 4). Library scoring complete (Section 5; binning/centroid internals open). RI calibration complete (Section 6; wizard auto-assignment open). Quantification ion selection complete (Section 7.1-7.3). Remaining quantification/statistics not analyzed.
