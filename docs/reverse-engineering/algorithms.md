# MetaboliteDetector 3.3 (Linux x86-64) - algorithm notes

Status: partial. Noise/baseline, peak detection and peak metrics are reconstructed. Deconvolution, library search, RI, quantification, MID, normalizers, statistics and PCA are NOT analyzed yet (see Section 7).

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
- Default values (Unknown): no default assignments were found in the GCMSSettings ctor (syms 0x3dfb98 / MCP 0x4dfb98 is empty) or in the static initializers. `GCMSSettingsObject` ctor (0x3f79de) only sets constructors for LS_MASS_FILTER and QUANT_EXCL_IONS. The values most likely come from a settings file on disk and/or GUI defaults. The settings storage is in static members (`_ZN4gcms12GCMSSettings<KEY>E`, syms 0xcc20a8..0xcc2170, plus 0xcc20e0 etc.), which are not initialized in the code I read. Next step: find `fromDisk` callers and any shipped default settings file.

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
  line through (minScanL, minL) and (minScanR, minR); abs = the lower of the two points
  repeat:
      for any i with (i*slope + icpt) - v[i] > 0.1: re-anchor the line through (abs, v[abs]) and (i, v[i])
  until no point violates
  residual r_i = v[i] - (slope*i + icpt) for all points, map<residual, index>
      (std::map insert: duplicate keys keep the first entry)
  take the first floor(N/2)+1 entries in ascending residual order
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
- Peaks created by DefaultPeakDetector get nf = getNoiseFactor() (Observed). The value used in the detector is set from analyzeBaseline's result `nf`. That linkage (whether detectPeaks writes stat.nf into this->nf) is Unknown.

## 3. Peak detection

### 3.1 DefaultPeakDetector<int,float>::detectPeaks (MCP 0x6ef020; syms 0x5ef020)
Setup (Observed):
```
if N < 5: return empty
BEGIN = AN_PEAK_THRESHOLD_BEGIN ; END = AN_PEAK_THRESHOLD_END
D  = SavitzkyGolayFilterFirstDerivative<float>(float(values))   # window 5 Inferred from size<5 check
S  = SavitzkyGolayFilter<int>(values)                         # window 5 Inferred
stat = analyzeBaseline(S, D, 15)                              # Section 2.1
```
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
- The exact handling of `i` after the tail branch (i = i-1 at one path) was recovered from assembly and may be off by one. Verify with a test case.
- Behavior when `i` equals the break index is Inferred.

### 3.2 DefaultPeakDetector<int,float>::createPeak (MCP 0x6efa1c; syms 0x5efa1c)
```
createPeak(ion, begin, end, max, values, D, extend, stat):
  if extend:
      while begin-1 valid and v[begin-1] < v[begin] (monotone toward baseline) and |v[begin] - stat.mean| >= 3*stat.sd: begin -= 1
      while end+1 valid  and v[end+1]  < v[end]   and |v[end]   - stat.mean| >= 3*stat.sd: end += 1
  apex = argmax over [begin, end] (first maximum wins, strict >)
  return Peak(ion, begin, end, apex, values[begin..end], nf = getNoiseFactor())
```
- The exact border conditions (direction check, and which side of stat is used) come from the decompile. Monotonicity check is Observed; exact comparison operator direction: Inferred.
- The "3 * sd" constant: Observed.

### 3.3 MatchedPeakDetector<int,float>::detectPeaks (MCP 0x6ee9ce; syms 0x5ee9ce)
```
W = peak_width (ctor arg); thr = ctor threshold (int)
f = tools::Filter::matched2ndDerivativeGaussianFilter<int>(values, 4*W)   # Observed call
regions: begin at i where f[i-1] < thr < f[i]; end at first j > i with f[j] <= thr
accept region if (end - begin) > 2
createPeak(ion, begin, end, f[begin..end])
```
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
- Observed (summary of decompile); the right-side condition should be verified with the decompile before porting.

## 4. Deconvolution / component detection

All addresses below are syms addresses (MCP = +0x100000). Ghidra cuts most of these functions at `std::endl` and `operator new` (it marks both no-return). The control flow after those calls was read from objdump.

### 4.0 Pipeline and entry points (Observed)
- GUI path: `metabolitedetector::ComponentDetector::detectComponents()` (0x4b3a0c) calls `GCMSAnalyzer::detectAllPeaks(reorder=true, scan_dist=10.0f, ion_dist=0.05f)` (call at 0x4b3a87; floats at rodata 0xa81010 = 10.0 and 0xa8100c = 0.05). Then it calls `deconvoluteChromatogram()` (0x5e907e).
- RI calibration wizard: `RICalibrationWizard::CompoundDetector::detect()` calls `detectAllPeaks(reorder=false, ...)` (0x3f4710), then `deconvoluteChromatogram()`.
- Batch quantification: `ResultSeries::scanIons` (0x638589..0x6385d7) builds an analyzer with a fresh `DefaultPeakDetector`, calls `setMinimumPeakNumberRequired(1)`, then `deconvoluteChromatogram(peaks, scan_count, start_scan, 2 * AN_DECONVOLUTION_WIDTH)` on a scan window.
- `MIN_PEAK_NUMBER` is field +0x0 of GCMSAnalyzer. The settings dialog sets it from `AN_MIN_PEAK_NUMBER` (call at 0x3c1177). The ctor default is Unknown.
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
f = tools::Filter::matched2ndDerivativeGaussianFilter<double>(v, W)   # filter internals: Section 7 gap 5
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
- AN_MIN_HEIGHT_BP is a fraction of the base peak of the deconvoluted spectrum. Its default value is Unknown (it is not in GCMSSettings::toStream, see Section 1).
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
- `tools::Filter::matched2ndDerivativeGaussianFilter<double>` kernel (sigma vs. W, kernel length, normalisation, edge handling). This is needed to reproduce the threshold T = 5.0 exactly.
- Default values of AN_DECONVOLUTION_WIDTH, AN_NUMBER_BINS, AN_MIN_PEAK_NUMBER, AN_MIN_HEIGHT_BP and the GCMSAnalyzer ctor default of MIN_PEAK_NUMBER.
- Meaning of `getDataType() == 2`, and how `Compound(rt, ints, ions, res)` uses `res`.
- `scan->vtbl[3]` (scan count) and `vtbl[4]/[5]` (ion count / EIC) identities are Inferred from usage.
- Whether the duplicated bin peaks from 4.5 are intended. A port should reproduce them first to match golden output.

## 5. Library search and RI - NOT ANALYZED
- RICalculator::calculateRI<int,float> at syms 0x3f98f4 (not decompiled beyond the entry).
- LibrarySearch, SpectrumFunctions: not decompiled.
- Settings in play: LS_RI_DIFF, LS_PURE_FACTOR, LS_IMPURE_FACTOR, LS_SCALE_FACTOR_ION/VALUE/NO_MATCH, LS_RESOLUTION, LS_USE_RI, LS_USE_SCALED_LIB, LS_USE_SPEC, LS_CUTOFF (names only).

## 6. Quantification, MID, normalizers, Statistics, PCA - NOT ANALYZED
- Settings in play: QUANT_ION_NUMBER, QUANT_MIN_DIST_IONS, QUANT_QUALITY_CUTOFF, QUANT_EXCL_IONS.
- tools::Statistics::collectedMean (syms 0x6c3a2c), collectedVariance (0x6c383a / 0x6c3640), polynomialRegression (0x6c3292): entry points only.

## 7. Open gaps (ranked by impact on reproducing results)
1. Deconvolution and component detection (Section 4). Without it no component table can be reproduced.
2. Library search scoring and RI (Section 5). Needed for identification output.
3. Default setting values. Storage format is known (Section 1); the default values are not in code I read. Need to find `fromDisk` callers and any shipped default settings file, or the GUI defaults.
4. DefaultPeakDetector details: exact SG filter coefficients and window (only 5 inferred), derivative definition, the C=0.2 tail rule (value read from rodata, not verified against output), the `i` update after a rejected peak (possible off-by-one), and the identity of vtbl[4]/[5].
5. Matched filter internals (matched2ndDerivativeGaussianFilter<int>), and whether MatchedPeakDetector should use raw or filtered values in createPeak.
6. Pooled statistics formulas (collectedMean/collectedVariance) and the meaning of the `n` field in BaselineStatistics.
7. Whether DefaultPeakDetector writes stat.nf into the noise factor used by Peak (Section 2.3).
8. Peak sharpness right-side condition (Section 3.8).
9. Quantification, MID, normalizers, statistics tests, PCA: not started.

## 8. Ranked next steps
- Decompile deconvoluteChromatogram and detectComponents (use objdump for truncated parts).
- Read LibrarySearch/SpectrumFunctions scoring and RICalculator::calculateRI.
- Locate `GCMSSettings::fromDisk` callers and any default settings file in the install.
- Verify the Section 3.1 loop against a real chromatogram output (golden data) before the Rust port.
