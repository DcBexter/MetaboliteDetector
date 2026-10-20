# MetaboliteDetector 3.3 legacy file formats (reverse-engineering notes)

Scope: formats written/read by the legacy Qt/C++ binary `metabolitedetector_app` (x86-64 ELF, PIE).
Method: Ghidra 12.1.4 decompilation via the rea-wsl MCP, DWARF from the unstripped build, nm symbols,
and byte-level checks on the sample library. No binary or library data is copied into the repository.

Conventions
- Every claim is tagged **Observed** (decompiled or byte-verified; the function and address are cited)
  or **Inferred** (derived from symbol names, local variable names, string literals, or symmetry).
- Address convention: Ghidra address = nm address + 0x100000 (PIE image base). Example:
  `GCMSSettings::toStream` nm `0x6d5c1c` = Ghidra `0x7d5c1c`. The addresses below are Ghidra addresses.
- Ghidra drops the length argument of many `std::ostream::write` calls and truncates some functions after
  calls to `std::string::c_str`, `operator new`, `istream::read` or `ostream::operator<<(endl)`, which Ghidra
  marks as no-return. Where a decompile is truncated this is stated and the claim is downgraded to Inferred.
- All integers are little-endian. "u32" = uint32, "i32" = int32, "f32"/"f64" = IEEE float/double.
- Strings written by `std::string` as `u32 length + raw bytes` (no terminator). Encoding is unverified.

---

## 1. Compound library (.lbr)

Container: gzip (`ogzstream`, Observed in `LibraryCompound::toDisk` 0x721070). Sample
`lib_aa.lbr` is gzip; `zcat` yields `lib_aa.raw` (152459 bytes, version 2, 27 compounds, fully consumed).

Top level (Observed, `LibraryCompound::toDisk` 0x721070): `ogzstream` open, then `u32 version = 5`
(writer), `u32 compound_count`, then `compound_count` x `LibraryCompound::toStream`.
The sample file has `version = 2` (produced by an older writer). The old reader
(`LibraryCompound::fromStream` 0x71fe5c) accepts versions 2 to 5; the gates are listed in 1.2.

Reference files: `tools/golden/make-lbr-v5.sh` makes version 5 files with the old program itself (reading
the sample and saving it, adding fragments with `addFragment`, an MSL export and import, and compounds from
a deconvoluted run) and dumps what the old reader gets back (`old.tsv`). They are never committed.

### 1.1 Top-level fields

| Offset | Type | Meaning | Tag |
|---|---|---|---|
| 0 | u32 | version (writer 5; sample 2) | Observed (toDisk 0x721070; sample and v5 files byte-checked) |
| 4 | u32 | compound count | Observed (toDisk 0x721070) |
| 8 | record | LibraryCompound records (see 1.2) | Observed |

### 1.2 LibraryCompound record (per compound)

Observed in `LibraryCompound::toStream` 0x71f832 and `fromStream` 0x71fe5c (disassembly), checked byte
by byte on the version 2 sample and on version 5 files written by the old program (`md-core/tests/lbr_sample.rs`,
2026-10-20: every field, peak and model value of 97 compounds matches `old.tsv`).

| # | Type | Meaning | Versions |
|---|---|---|---|
| 1 | Spectrum | reference spectrum (see 1.4) | all |
| 2 | Compound | compound record (section 3.2): n_peaks, peaks, model, name ... check | all |
| 3 | u32 count + (f32 mz, i32 intensity)[count] | `comb_spec` (combined spectrum, `map<float,int>`) | all |
| 4 | u32 count + (f32, u32)[count] | `ion_freq` (ion frequency map) | all |
| 5 | u32 | `spec_count` | all |
| 6 | u32 count + f64[count] | `ris` (retention index list) | all |
| 7 | u32 count + 40-byte peak records (1.3) | old peak list, given to the compound with `setPeaks` | 2 and 3 only |
| 8 | u32 count + string[count] | `sources` | all |
| 9 | u32 count + (u32 begin, u32 end, string formula)[count] | `formulas` (`pair<ulong,ulong>` written as 2 x u32); the reader trims each formula (`trimFormulas` 0x721876, boost `trim`) | 3 and later |
| 10 | u32 count + (u32 begin, u32 end, u32 size)[count] | `mid_sizes` | 5 only |

Both maps are sorted by key: `addFragment(73,147," C3H9Si ",2)` then `(5,6,"N",0)` is written as
`(5,6,"N"), (73,147,"C3H9Si")` and mid sizes `(5,6,0), (73,147,2)` (Observed, `fragments_v5.lbr`).

Versions 2 and 3 compounds with `n_peaks > 0` store `Peak<int,float>` records whose layout is not decoded;
the rewrite rejects them. The sample and all v5 files have none of those (Observed).

### 1.3 Peak record of versions 2 and 3 (40 bytes)

`PeakEssential<float>` in its pre-version-4 form (Observed: `PeakEssential<float>::fromStream` 0x6fa23c reads
these fields and stops for versions below 4; the old program turns them into the v5 peaks of 1.6).

| Offset in record | Type | Meaning | Value in sample |
|---|---|---|---|
| 0 | f32 | ion (m/z) | m/z |
| 4 | f64 | sharpness | 1.0 |
| 12 | f64 | integral | intensity at that m/z |
| 20 | f64 | quality | 0.0 |
| 28 | f64 | signal to noise | -1.0 |
| 36 | u32 | count | 1 |

Saved as version 5, the remaining fields get first/last scan 0, model parameter 1, max scan 0 and
interpolated max scan 0 (Observed, `lib_aa_v5.lbr`). The rewrite keeps these records as raw `blocks`.

### 1.4 Spectrum (used in 1.2 #1)

Observed: `Spectrum<int,float>::toStream` 0x7027a2 (decompile) and byte check in `parse_lbr.py`.

| Offset | Type | Meaning |
|---|---|---|
| 0 | u32 | scan number |
| 4 | f64 | rt in ms (Observed: the MSL export divides it by 60000 for minutes, 6.2) |
| 12 | f64 | ri |
| 20 | f32 | bin_width |
| 24 | u32 | n_masses |
| 28 | f32[n_masses] | masses |
| ... | u32 | n_values |
| ... | i32[n_values] | values (intensities) |
| ... | u32 | check = (scan + n_masses) mod 2^32 (Observed in sample) |

### 1.5 Validated parsers

`docs/reverse-engineering/tools/parse_lbr.py` (`python3 -I parse_lbr.py < file.lbr`) decodes version 2 only.
The Rust reader `md_core::library::parse_lbr` reads versions 2 to 5 and returns the same values for every
field of all 27 sample compounds (`tools/dump_lbr.py` against `md-core/tests/lbr_sample.rs`, 2026-10-19) and
of all v5 reference files (1, 2026-10-20). It rejects version 2 and 3 compounds with peak records (1.2).

### 1.6 Peaks and model peak (versions 4 and 5)

Written inside the Compound record (3.2) after `n_peaks`. Observed in `PeakEssential<float>::toStream`
0x6f994c, `PeakEssential<short>::toStream` 0x6facfe, `Peak<double,short>::toStream` 0x6f66e2 and their
`fromStream` (0x6fa612, 0x6fb9c7, 0x6f680b); member offsets from the disassembly, values checked on the v5 files.

`PeakEssential<T>` (68 bytes for T = f32, 66 for T = i16), in file order:

| # | Type | Meaning | Member |
|---|---|---|---|
| 1 | T | ion | +0x08 |
| 2 | f64 | sharpness | +0x10 |
| 3 | f64 | integral | +0x20 |
| 4 | f64 | quality | +0x18 |
| 5 | f64 | signal to noise | +0x28 |
| 6 | u32 | count (`size_t` in memory) | +0x30 |
| 7 | u32 | first scan | +0x38 |
| 8 | u32 | last scan | +0x3c |
| 9 | f64 | model parameter | +0x50 |
| 10 | u32 | max scan | +0x40 |
| 11 | f64 | interpolated max scan | +0x48 |

- Peaks: `n_peaks` x `PeakEssential<float>`, sorted by ion (the compound keeps `map<float, PeakEssential*>`).
- Model (always present from version 4): `Peak<double,short>` = `PeakEssential<short>`, `f64 noise_factor`,
  `u32 n`, `f64[n]` values, `u32 check = first + last + max + n` (mod 2^32; the reader throws on mismatch).
  Without a model the old program writes ion 0, sharpness -1, integral -1, quality 0, s/n -1, count 0,
  scans 0, model parameter 1, noise factor 0.5 and one value 0 (Observed in `lib_aa_v5.lbr`).

### Unknowns (section 1)
- `Peak<int,float>` layout of version 2 and 3 compounds with `n_peaks > 0` (no such file seen).
- String encoding (latin-1 assumed by the parsers; the real encoding is unverified).

---

## 2. Chromatogram data (.bin + .idx)

Symbols: `GCMSDiskScanWriter<int,float>` (writer), `GCMSDiskScan<int,float>` (reader, wraps a
`DataProvider<int,float>`), `DataProvider<int,float>::readIndexFile` / `writeIndexFile`.

Sample files: `docs/reverse-engineering/tools/diskscan/make-cache.sh` runs the old binary's own
`CdfImporter(cdf).importData(stem)` on the 15 example runs (`LD_PRELOAD` driver). The output is
deterministic (two runs give byte-identical files). Everything marked "byte-verified" below was checked
with a parser on all 15 pairs: every index is consumed exactly and every record header matches its
index entry. The reader is `md_core::diskscan`; the files are never committed.

File naming (Observed)
- `importData(stem)` writes `<stem>.bin` (data) and `<stem>.idx` (index) (byte-verified).
- `DataProvider::readIndexFile` 0x75bf30 appends ".idx" to `f_name`; `GCMSDiskScanWriter::writeIndexFile`
  0x6e7d00 appends ".idx" to `cache_file` and opens `out|trunc|binary`.
- `CdfImporter::importData` 0x7cc4cc creates the writer with `DataType = 2` (Observed: `mov $0x2,%edx`
  before the `GCMSDiskScanWriter(char const*, DataType)` call at 0x7cc512).
- Every scan is centroided before it is written (`CdfImporter::scanAdd`, algorithms.md §4.1). The stored masses are nominal: on a synthetic run with exact masses all 4400 scans hold the centroids rounded to whole numbers.

### 2.1 Scan record (data file)

Written by `GCMSDiskScanWriter::addScan(uint scan, double time, vector<float>& masses,
vector<int>& values, streampos pos, double tic)` 0x6e6888 (Observed; byte-verified).
The vector overload 0x6e7072 calls it with `pos = 0` and `tic = -1.0` (Observed), so that path writes the
masses. The other overloads (0x6e6fd6, 0x6e65c8) were not decompiled.

| Offset in record | Type | Meaning | Tag |
|---|---|---|---|
| 0 | i32 | `size_all`: total record length in bytes including this field. Masses written: 36 + 8n. Masses shared: 36 + 4n | Observed, byte-verified |
| 4 | i64 | `masses_ptr`: absolute file offset of the mass array (record start + 36 + 4n when masses written, else the given pos) | Observed, byte-verified |
| 12 | f64 | `tic`: total ion current; `CdfImporter` passes the `.cdf`'s `total_intensity`, which is not the sum of the spectrum | Observed, byte-verified |
| 20 | u32 | scan number (0-based) | Observed, byte-verified |
| 24 | f64 | retention time in ms. `CdfImporter` stores `(double)(float)seconds * 1000`, so times carry f32 precision | Observed, byte-verified |
| 32 | u32 | n = number of values (ions) in this scan | Observed, byte-verified |
| 36 | i32[n] | intensities (`s_values`) | Observed, byte-verified |
| 36 + 4n | f32[n] | masses (`s_masses`), present only when masses are written | Observed, byte-verified |

Header of the data file: none. The first record starts at offset 0 and records follow each other
(byte-verified).

Side effects of `addScan` (Observed, 0x6e6888), kept in memory and written to the index:
- `recptr[scan]` = file offset of the record start (i64)
- `massptr[scan]` = `masses_ptr` (i64)
- `scan2tic[scan]`, `scan2rt[scan]` (f64), `scan2numval[scan]` (int, n)
- `ions`: set of all distinct mass values (float)
- `scan2stat[scan]`: BaselineStatistics (only for COMPRESSED_PROFILE and CENTROID types)

### 2.2 Compression

- Type `COMPRESSED_PROFILE` (enum value not recovered): `analyzeBaseline(values, 100)` produces a
  BaselineStatistics stored per scan, then `compressSpectrum` runs (Observed: 0x6e6888 -> 0x6e758c).
- `compressSpectrum` (Observed, 0x6e758c): if fewer than 10000 masses, prints "Only N values..." to
  `std::cout` (debug). Calls `SpectrumFunctions::detectPeaks<int>(peaks, s_values, thresh_beg, thresh_end,
  base_dev, true, stat)`, then keeps only the intensities inside detected peaks and the corresponding masses
  (taken from `s_masses[begin + i]`), replacing both vectors. Thresholds come from
  `GCMSDiskScanWriter::setCompressionParameter(double, double, double)` 0x6e7164 (values not recovered).
  The record layout stays the same, so a reader needs no special case.
- `CENTROID` type (Observed, 0x6e6888): stores a zeroed BaselineStatistics and writes the data unchanged.
- Other profile types: data written unchanged (Observed, 0x6e6888: no branch).
- The CDF import (`DataType = 2`) writes the spectra unchanged with zeroed statistics (byte-verified:
  the spectra equal the `.cdf`). Whether 2 is CENTROID is Inferred.

### 2.3 Index file (.idx)

Reader: `DataProvider::readIndexFile` 0x75bf30 (Observed from the disassembly: the sizes of the
`istream::read` calls, the checks and the map inserts; Ghidra's decompile is truncated). Byte-verified on
all 15 sample files. Writer: `GCMSDiskScanWriter::writeIndexFile` 0x6e7d00 (not disassembled; the
sample files match the reader).

| Item | Type | Meaning |
|---|---|---|
| magic | i64 | `0x4b494104`; any other value throws "Error opening index file!" |
| size_int, size_double, size_long_long | u32 x3 | must be 4, 8, 8, else throws "Number of used bytes for primitive data types does not match." |
| scan count | u32 | number of per-scan entries below |
| mass size | u32 | byte size of one mass (4 = float); also the stride of the ion list |
| value size | u32 | byte size of one intensity (4 = int) |
| feature count | u32 | then per feature: u32 length + key bytes, u32 length + value bytes. The CDF import writes the 11 ANDI attributes it looks for (`acq_time`, `experiment_type`, `number_of_scans`, ...); missing ones as " " |
| TOF polynomial | u32 count + f64[count] | `FilteredDataProvider` applies it to the masses; the CDF import writes count 0 |
| global baseline statistics | 4 x f64 + i64 | zero in the samples; order follows `BaselineStatistics(double, double, double, double, long long)` |
| RI calibration | `RICalculator::fromStream` (section 4.1) | the run's calibration from the RI wizard; samples: 0 entries, drt 0.0, source "not calibrated" |
| ion count | u32 | then `ion count * mass size` bytes of distinct masses; the reader skips them with `seekg` |
| per-scan entries | scan count x 84 bytes | see below |

Per-scan entry (84 bytes):

| Offset | Type | Meaning |
|---|---|---|
| 0 | u32 | scan number |
| 4 | i64 | file offset of the record in the `.bin` |
| 12 | i64 | `masses_ptr` (same as in the record) |
| 20 | f64 | tic |
| 28 | f64 | retention time in ms |
| 36 | u32 | n (values in the scan) |
| 40 | u32 | DataType, copied into the provider for each entry (2 in the samples) |
| 44 | 4 x f64 + i64 | BaselineStatistics of the scan (zero in the samples) |

The reader fills `rt2scan`, `scan2rt`, `recptr`, `massptr`, `scan2numval` and `scan2stat` from the
entries and finally throws "Error reading data." if the stream failed.

### 2.4 Reading

`GCMSDiskScan::getSpectrum(unsigned int)` 0x6e2638 (truncated): calls `DataProvider::moveToScan`,
`getTime`, `getIonCount`, `getSignals(0, n-1)`, `getMasses(0, n-1)`. With the index above, a reader
seeks to the record offset for the intensities and to `masses_ptr` for the masses. `md_core::diskscan`
does that and checks each record header against its index entry.

### Unknowns (section 2)
- `DataType` enum names for the values. 2 is what `CdfImporter` writes (nominal masses); `FastFlightImporter` writes 1 or 0 (exact masses), algorithms.md §4.1.
- `setCompressionParameter` default thresholds and the meaning of the BaselineStatistics fields.
- Whether `writeIndexFile` can write another magic or other type sizes (only the reader was disassembled).
- None for the TOF polynomial: `FilteredDataProvider` maps `m/z = (f32) poly(bin as f64)`, which `md_core::diskscan` applies (section 4.3).

---

## 3. Detected compounds (.cmp)

Symbols: `Compound<int,float>`, `Spectrum<int,float>`, `PeakEssential<float>`, `Peak<double,short>`,
`PeakEssential<short>`, `Peak<int,float>`.

Verified on 15 reference files the old 3.3 binary wrote itself (`tools/golden/make-reference.sh` with
`MD_RESULTS`, never committed): every field below decodes, every value matches the old program's own
getters (`raw.tsv`, 9 digits) and `md_core::cmp` writes the payload back byte for byte. The old program
reads files written by `md_core::cmp` (`make-reference.sh --reread`): rewriting them gives the same payload
except `time_stamp`, which it also changes when it rereads its own files (Observed). A component with an empty
name (md-core writes none: the rewrite has no name for it) gets a 15-digit number as name when the old program
reads it, and the check value (§3.2, item 13) follows (Observed on two example runs, 34 components).

### 3.1 File container

- `Compound<int,float>::toDisk(vector<Compound*>, file)` 0x711332 (Observed): `ogzstream` (gzip) around
  `compoundsToStream` 0x7111d4: `u32 version = 4`, `u32 count`, then `Compound::toStream` per compound.
- No file-level checksum; each spectrum, model peak and compound ends with a u32 check value.

### 3.2 Compound record (`Compound<int,float>::toStream` 0x70fcba)

Order from the disassembly (Observed) and byte-verified (Verified):

| # | Type | Meaning |
|---|---|---|
| 1 | Spectrum record | `Spectrum<int,float>::toStream` as in `.lbr` (§1.2): u32 scan, f64 rt (ms), f64 ri (-1 without calibration), f32 bin width, u32 n + f32 masses, u32 n + i32 values, u32 check = scan + n |
| 2 | u32 n + PeakEssential<float>[n] | ion peaks behind the spectrum (`std::map<float, PeakEssential*>`, ascending ion), 68 bytes each (3.3) |
| 3 | Peak<double,short> | model peak (`unique_ptr`, written unconditionally), 3.3 |
| 4 | u32 length + bytes | name (no terminator, no encoding) |
| 5 | i32 n + f64[n] | fwhm list; n = -1 for null |
| 6 | i32 n + f64[n] | integrals list; n = -1 for null |
| 7 | u32 n + f32[n] | quant_ions (`std::set<float>`) |
| 8 | u32 n + f32[n] | auto_quant_ions (`std::set<float>`) |
| 9 | u32 n + (string, string)[n] | features (`std::map<string,string>`) |
| 10 | f64 | drt |
| 11 | i64 | time_stamp (`time_t`) |
| 12 | i64 | id |
| 13 | u32 | check = len(quant_ions) + len(features) + len(name bytes) + n peaks (mod 2^32) |

### 3.3 Peak records

`PeakEssential<float>::toStream` 0x6f994c (Observed write sizes; field meaning Verified against the getters):

| Offset | Type | Meaning |
|---|---|---|
| 0 | f32 | ion (m/z) |
| 4 | f64 | sharpness |
| 12 | f64 | integral |
| 20 | f64 | quality |
| 28 | f64 | signal to noise |
| 36 | u32 | count (meaning not recovered; 1 in every sample) |
| 40 | u32 | first scan |
| 44 | u32 | last scan |
| 48 | f64 | model parameter (share of the model peak in the trace) |
| 56 | u32 | apex scan |
| 60 | f64 | interpolated (precise) apex scan |

`PeakEssential<short>::toStream` 0x6facfe: the same fields, but the ion is an i16 (66 bytes). The model
peak's ion is 1 in the samples.

`Peak<double,short>::toStream` 0x6f66e2: `PeakEssential<short>`, f64 noise factor (`setNoiseFactor`
writes the same member, offset 0x58), u32 n + f64[n] model values (one per scan from first to last
scan), u32 check = n + first scan + last scan + apex scan.

These are the peak records of version 4 and 5 libraries too (§1.6); `md-core` reads both with one
reader (`library::read_compound_record`).

`Peak<int,float>::toStream` 0x6f3036 (not used by 3.3 `.cmp` files): `PeakEssential<float>`, f64,
u32 n + i32[n], u32 check = n + first + last + apex scan.

### 3.4 Reader and versions

`Compound<int,float>::fromStream(istream&, Compound&, version)` 0x712ece (Observed branches): after the
spectrum it reads the peak count; version > 3 reads `PeakEssential<float>` records and the
`Peak<double,short>` model peak, version <= 3 reads `Peak<int,float>::fromStream(version)` records and no
model peak. `md_core::cmp` reads version 4 only and reports other versions as unsupported.

### Unknowns (section 3)
- Layout of `.cmp` versions 1 to 3 (older MetaboliteDetector releases); no sample file.
- Meaning of the `count` field and of the model peak's ion code.
- Where the reader sets `time_stamp` to the load time (seen only in its effect).

---

## 4. Calibration files

### 4.1 RI calibration (.calib): `gcms::RICalculator`

Members: `rt2ri` (`std::map<double,double>`), `drt` (f64) and `source` (`std::string`, `getSource`).
Verified with a calibration the harness built through `RICalculator(name)` and `addCalibrationPair(rt, ri)`
and saved with `toDisk` and `toStream` (`tools/golden/make-reference.sh`, `MD_RESULTS`): the values are
stored as given, so the unit is the caller's. The GUI passes spectrum retention times, which are ms.

Binary form, `RICalculator::toStream` 0x7d2fd2 (Verified):

| Offset | Type | Meaning |
|---|---|---|
| 0 | u32 | entry count (`rt2ri.size()`) |
| 4 | f64 rt (ms), f64 ri | repeated count times, ascending rt (std::map order) |
| 4 + 16n | f64 | drt (ms; Inferred: a difference of spectrum RTs, algorithms.md §6.5; every sample has 0) |
| ... | u32 + bytes | source (no terminator); "not calibrated" in uncalibrated runs |

The binary form is not a file of its own: it sits in the run's `.idx` (§2.3). The RI wizard and the
calibration setter store a run's calibration there (`CalibrationPage::setDataToChromatogram`,
`RISetter::recalibrateChromatograms`, `CalibrationSetter::on_b_start_ri_clicked`: `setRICalculator`, then
`writeIndexFile`; Observed calls). `md_core::diskscan::read_calibration` reads it.

Text form (the `.calib` file), `RICalculator::toDisk(file)` 0x7d28f0 (Verified): one line per entry,
`<rt>` TAB `<ri>` and a newline, default ostream precision (6 significant digits, so `612345.678901`
becomes `612346` and `1500000` becomes `1.5e+06`). drt and source are not stored. The text form is lossy.
`fromDisk` 0x7d22ac splits lines on `\s+` (Observed) and reads the harness file back to the same text.

### 4.2 TOF calibration (.mzcalib)

Writer: `TOFCalibrationWindow::on_act_save_triggered` 0x6bc1cc (Observed; decompile complete).
Filter string "Calibration files (*.mzcalib)" and "TOF calibration (*.mzcalib);;All files (*.*)" in the binary.

- Opens the file with `QFile` in `WriteOnly | Truncate | Text` mode (CRLF line endings on Windows are
  Inferred from Qt text mode). Appends `.mzcalib` if the name does not already end with it (Inferred
  from `endsWith` and `remove(size - 6)`).
- Text lines via `QTextStream` (Observed, in order):

| Line | Format | Meaning |
|---|---|---|
| 1 | `Polynomial order:` TAB `<order>` | order = parameter.size() - 1 |
| 2..(1+N) | `Parameter:` TAB `<i>` TAB `<value>` | polynomial coefficient i, N = parameter.size() |
| next | `Calibration table` | fixed header |
| next..end | `<key>` TAB `<value>` | one line per `CalibrationPairWidget::getCalibrationPairs()` entry; key = `float` (m/z), value = `double` |

- Precision of the numbers: `setRealNumberPrecision(18)` (Observed, disassembly), i.e. C `%.18g`.
- Encoding: `QTextStream` default (Unknown; ASCII content only in the strings).

Reader: `MainWindow::readTOFCalibData(QString file, QWidget*)` 0x4c8aba (Observed, disassembly). A
whitespace-token `ifstream`: skips two tokens (`Polynomial order:`), reads `order` as int, then for `i` in `0..=order`
skips the label, reads the index and the double value. The calibration table is not read. Ported in
`md_core::tof::TofCalibration::parse` (an order outside 0..=64 or a missing/invalid value is an error).

Regression: `TOFCalibrationCalculator::regress` 0x6db8a2 needs at least `degree + 1` pairs (error text
"ERROR: Not enough values for polynomial regression."), then GSL `multifit_linear` (SVD) with
`X[i][j] = bin_i^j`, `y = m/z_i` (Observed). Ported as `TofCalibration::fit` (nalgebra SVD).

The window of the rewrite (`crates/md-app/src/tof_calibration.rs`, #148) sums a run per TOF bin
(`fastflight::bin_spectrum`), takes the TOF bin of a clicked peak as the intensity-weighted centre of the points
at or above half the local maximum (`tof::peak_bin`), asks for the reference m/z by hand, fits with
`TofCalibration::fit` and saves with `TofCalibration::to_text`. The old window's automatic calibration is not
ported (below), so a calibration needs reference ions chosen by hand.

`TOFCalibrationWindow::saveSettings()` 0x6c077e is **not** the .mzcalib writer. It stores only the window
size in `QSettings("TU-Braunschweig", "TOFCalib")`, group "GUI", key "size" (Observed).

### Unknowns (section 4)
- Encoding of the text files (ASCII numbers only in `.calib`).
- Automatic TOF calibration (`automaticCalibration`, `iterativeCalibration`, `searchBestMatch`): not decoded.

### 4.3 FastFlight raw data (MsData.FFT + MsData.FFC): `gcms::FastFlightImporter`

Observed (disassembly of `initialization`, `readSpectrum<float>`, `importData`; Ghidra address = nm + 0x100000).
No real FastFlight file was available: the layout is checked against a synthetic run only, read by the old
importer itself (`docs/reverse-engineering/tools/fastflight/`).

The importer takes a directory (trailing `/`; the names `MsData.FFT` and `MsData.FFC` are appended).

**FFC (index).** `i32 N`, then N entries of 36 bytes: 4 bytes skipped, `i32 n`, 4 skipped, `u64 tic`,
`u64 fft_offset`, `f64 factor`. Time of the scan = `factor * n` (factor 0 means 2.5); the unit is ms (Inferred).
A negative `n` or offset is "Error reading data!". Scans are keyed by time (the first of equal times stays).

**FFT (data).** No file header; a spectrum starts at the offset from the FFC. All words are 24 bit, little endian.
A frame is `0xFFFFFF`, a header word, then data words (else "Expected new frame header!"):
- header: count = `h & 0x3FFF` (words in the frame, including the header words); bit 14 = time block, bit 15 = final
  frame (tested only if bit 14 is clear). The first `skip` words of the count are header: 2, 3 with a time word
  (`& 0x1FFFFF`, times 2.5; clears the point lists), 4 for the final frame.
- data words: `0xFF8000` is followed by a base offset; `(w & 0xFF0000) == 0xFF0000` sets the position
  `w & 0xFFFF`; any other word is an intensity at TOF bin `base + position`, then the position increases by 1.
- after the final frame: `w1 + (w2 << 24)` as i32 is the TIC. The smallest intensity is subtracted from all.

**importData(path, dropFirst, b).** For each scan in time order: TIC of the FFT trailer must equal the FFC one
("TIC mismatch! TICs in the file .FFT and .FFC do not match."); TIC 0 is skipped; `dropFirst` (the batch importer
passes 1) removes the first point. With `b` false the scan goes through peak detection and compression; with `b`
true only scans with more than 50 000 points are written, as they are. The cache gets the raw TOF bins as
"masses" and the polynomial in the `.idx` (section 2.3); the reader maps `m/z = (f32) poly(bin as f64)`.
Ported in `md_core::fastflight`; the compressing path is not.
- Unknown: whether a time block also resets the minimum intensity; the port keeps it (a time block comes first in real spectra).
- The Rust reader also keeps small scans uncompressed; neither old path does that. It refuses an index whose scans share FFT data (more points than FFT words).
- The old writer's `b = true` cache repeats points (runs of 1 to 5) and loses about 0.5 % of them (found with the
  synthetic run); the Rust test compares as a sequence, not point by point.

---

## 5. GCMSSettings (SF_*, AN_*, LS_*, SQL_*, QUANT_*)

### 5.1 Default values

Static members (`gcms::GCMSSettings::*`) are initialized from `.data` (Observed: nm symbol addresses and
sizes; values read from the ELF `.data` with file offset = vaddr - 0x1000). The constructors
`GCMSSettings::GCMSSettings()` 0x4dfb98 and `GCMSSettingsObject::GCMSSettingsObject()` 0x4f79de assign no
scalar defaults (Observed). The static-init function 0x7d73ee only calls `std::ios_base::Init` (Observed).

Types: from symbol size and bit pattern (Inferred). The bool fields are 1 byte. The 4-byte fields are
either i32 or float; the bit pattern decides for LS_RESOLUTION (0x3e99999a = 0.3f).

| Field | Size | Type | Default | Tag |
|---|---|---|---|---|
| SF_PEAK_THRESHOLD_BEGIN | 8 | f64 | 10 | Observed bytes |
| SF_PEAK_THRESHOLD_END | 8 | f64 | -10 | Observed bytes |
| SF_MAX_BASE_DEV | 8 | f64 | 20 | Observed bytes |
| SF_FWHM | 8 | f64 | 0.5 | Observed bytes |
| SF_CENTROID_THRESHOLD_BEGIN | 8 | f64 | 10 | Observed bytes |
| SF_CENTROID_THRESHOLD_END | 8 | f64 | -5 | Observed bytes |
| SF_CENTROID_MAX_BASE_DEV | 8 | f64 | 30 | Observed bytes |
| AN_PEAK_THRESHOLD_BEGIN | 8 | f64 | 10 | Observed bytes |
| AN_PEAK_THRESHOLD_END | 8 | f64 | -5 | Observed bytes |
| AN_MAX_BASE_DEV | 8 | f64 | 20 | Observed bytes |
| AN_MIN_PEAK_HEIGHT | 8 | f64 | 20 | Observed bytes |
| AN_DECONVOLUTION_WIDTH | 8 | f64 | 1 | Observed bytes |
| AN_NUMBER_BINS | 4 | i32 | 10 | Observed bytes |
| AN_MIN_PEAK_NUMBER | 4 | i32 | 10 | Observed bytes |
| AN_MODELED_VALUES_FOR_INTEGRATION | 1 | bool | true | Observed bytes |
| AN_MIN_HEIGHT_BP | 8 | f64 | no static value (.bss zero-init) | Observed (.bss) |
| LS_RI_DIFF | 8 | f64 | 20 | Observed bytes |
| LS_PURE_FACTOR | 8 | f64 | 0.5 | Observed bytes |
| LS_IMPURE_FACTOR | 8 | f64 | 0.5 | Observed bytes |
| LS_CUTOFF | 8 | f64 | 0.8 | Observed bytes |
| LS_SCALE_FACTOR_ION | 8 | f64 | 2 | Observed bytes |
| LS_SCALE_FACTOR_VALUE | 8 | f64 | 0.6 | Observed bytes |
| LS_SCALE_FACTOR_NO_MATCH | 8 | f64 | 1 | Observed bytes |
| LS_RESOLUTION | 4 | f32 | 0.3 | Observed bytes (bit pattern) |
| LS_USE_RI | 1 | bool | true | Observed bytes |
| LS_USE_SCALED_LIB | 1 | bool | true | Observed bytes |
| LS_USE_SPEC | 1 | bool | true | Observed bytes |
| LS_MASS_FILTER | 0x30 | std::set<float> | empty (ctor-built) | Observed (.bss) |
| LS_IDENTIFICATION_ALGORITHM | 4 | enum | no static value (.bss zero-init) | Observed (.bss) |
| SQL_MAX_NUMBER_FRAGMENTS | 4 | i32 | 5 | Observed bytes |
| SQL_MINIMAL_NUMBER_EQUAL_FRAGMENTS | 4 | i32 | 3 | Observed bytes |
| QUANT_ION_NUMBER | 4 | i32 | 3 | Observed bytes |
| QUANT_MIN_DIST_IONS | 8 | f64 | 0.5 | Observed bytes |
| QUANT_QUALITY_CUTOFF | 8 | f64 | 1 | Observed bytes |
| QUANT_EXCL_IONS | 0x30 | std::set<pair<float,float>> | empty (ctor-built) | Observed (.bss) |

Note: `std::set` contents are built at runtime and are empty in the static image.

### 5.2 Serialization: GCMSSettings::toStream 0x7d5c1c (Observed)

Writes the static members in the order below. Offsets assume the default empty QUANT_EXCL_IONS.
Every `ostream::write` in the decompile matches one row; sizes from the symbol sizes above.

| Offset | Type | Field |
|---|---|---|
| 0 | f64 | SF_PEAK_THRESHOLD_BEGIN |
| 8 | f64 | SF_PEAK_THRESHOLD_END |
| 16 | f64 | SF_MAX_BASE_DEV |
| 24 | f64 | SF_FWHM |
| 32 | f64 | SF_CENTROID_THRESHOLD_BEGIN |
| 40 | f64 | SF_CENTROID_THRESHOLD_END |
| 48 | f64 | SF_CENTROID_MAX_BASE_DEV |
| 56 | f64 | AN_PEAK_THRESHOLD_BEGIN |
| 64 | f64 | AN_PEAK_THRESHOLD_END |
| 72 | f64 | AN_MAX_BASE_DEV |
| 80 | f64 | AN_MIN_PEAK_HEIGHT |
| 88 | f64 | AN_DECONVOLUTION_WIDTH |
| 96 | i32 | AN_NUMBER_BINS |
| 100 | i32 | AN_MIN_PEAK_NUMBER |
| 104 | f64 | LS_RI_DIFF |
| 112 | f64 | LS_PURE_FACTOR |
| 120 | f64 | LS_IMPURE_FACTOR |
| 128 | f64 | LS_SCALE_FACTOR_ION |
| 136 | f64 | LS_SCALE_FACTOR_VALUE |
| 144 | f64 | LS_SCALE_FACTOR_NO_MATCH |
| 152 | f32 | LS_RESOLUTION |
| 156 | u8 | LS_USE_RI |
| 157 | u8 | LS_USE_SCALED_LIB |
| 158 | u8 | LS_USE_SPEC |
| 159 | i32 | QUANT_ION_NUMBER |
| 163 | f64 | QUANT_MIN_DIST_IONS |
| 171 | f64 | QUANT_QUALITY_CUTOFF |
| 179 | u32 | n = QUANT_EXCL_IONS.size() |
| 183 | (f32, f32)[n] | QUANT_EXCL_IONS entries (second value written through a double variable; width 4 Inferred from the pair type) |
| 183 + 8n | f64 | checksum = n + SF_PEAK_THRESHOLD_BEGIN + LS_PURE_FACTOR + AN_PEAK_THRESHOLD_BEGIN (Observed formula) |

Total size: 191 + 8n bytes. For defaults with no excluded ions: 191 bytes. The checksum is 20.5 for defaults.

Not serialized by `GCMSSettings::toStream` (Observed by absence in the write sequence):
AN_MODELED_VALUES_FOR_INTEGRATION, AN_MIN_HEIGHT_BP, LS_CUTOFF, LS_MASS_FILTER, LS_IDENTIFICATION_ALGORITHM,
SQL_MAX_NUMBER_FRAGMENTS, SQL_MINIMAL_NUMBER_EQUAL_FRAGMENTS.

`SQL_MAX_NUMBER_FRAGMENTS` and `SQL_MINIMAL_NUMBER_EQUAL_FRAGMENTS` (Observed, disassembly) are set in the
settings dialog (`SettingsDialog::updateGUI` / `applySettings`) and stored by the app's QSettings under the
same key names (`MainWindow::saveSettings` / `readSettings`), not in this file. `LS_IDENTIFICATION_ALGORITHM` (Observed, disassembly): only two loads exist, in `getGCMSSettingsObject`
(0x6d683a) and in `SQLLibrary::getLibraryHits` (0x6e30cc); no instruction stores to it or takes its
address, so it keeps its zero-initialised `.bss` value 0 in every session. 0 is the path ported and verified
in §6.3 (the two spectrum scores, vtable +0x18 and +0x10). 1 scores fragments (`getSpectrumFragments` of
the query with bin width 2, then one virtual score call), 2 a dot product (`scaleCompound_dot` on both
spectra, then mass overlap counts); both are unreachable in the shipped app and are not ported.

Note: `LS_CUTOFF` (default 0.8) and the SQL fields are not persisted by this path. The brief asks whether
they are stored elsewhere; this pass did not check that (Unknown).

### 5.3 GCMSSettingsObject (instance copy) 0x7d6a1a (Observed)

The same field sequence as 5.2 written from instance members. One exception: the write for
AN_MAX_BASE_DEV uses the static address (`0xdc20f0` in the decompile), not the instance member. This is
probably an artifact of how the source compares or the global was used in the original code. Treat it as
a possible bug in the legacy writer (Inferred from the decompile). The checksum uses instance members.

### 5.4 File forms

- `GCMSSettings::toDisk(filename)` 0x7d64c6 (Observed): `ofstream(filename)`, then `toStream`, then close.
  The file has **no header, magic or version**. It is the raw 5.2 stream.
- `GCMSSettings::fromDisk(filename)` 0x7d657a (Observed): creates a temporary `GCMSSettings s`, opens
  `ifstream`, and if the stream is good calls `fromStream` (static), then closes. The temporary `s` is not
  used for loading (the static members are loaded). Whether truncated files are rejected is not verified
  (fromStream decompile is truncated).
- `GCMSSettings::fromStream(istream&)` 0x7d6052: decompile truncated at `istream::read` (Observed).
  The read order is Inferred to mirror 5.2 (locals `tmp`, `size`, `begin`, `end`, `i`, `check` match the
  write sequence). The checksum comparison is not visible.
- `GCMSSettingsObject::fromStream` / `fromDisk`: not decoded (Gap).

### Unknowns (section 5)
- Whether `fromStream` validates the checksum and whether it accepts files without the trailing set.
- Persistence of LS_CUTOFF and AN_MODELED_VALUES_FOR_INTEGRATION (not in this stream; SQL_* are QSettings).
- Exact type of LS_RESOLUTION (f32 inferred); LS_IDENTIFICATION_ALGORITHM is compared as a 32-bit int.
- Whether the AN_MAX_BASE_DEV static write in GCMSSettingsObject is intentional.

---

## 6. MSL (AMDIS) and NIST MSP import/export

Status: MSL partially decoded (the import function is decompiled but truncated; most of 6.1 is based on string
literals and local names, Inferred). NIST text decoded with a black-box driver (6.3).

### 6.1 MSL import: `LibrarySearch<int,float>::importMSL(char*)` 0x72b098 (Observed decompile, truncated)

Locals (Observed): `exp`, `exp_quant`, `exp_feature` (three `std::regex` objects), `name`, `ri`, `rt`,
`resolution` (initialized to "1.0"), `quant_ions` (set<float>), `features` (map<string,string>), `source`,
`line`, `u_line`, `value_1`, `match`, `ion`, `f_res`, `d_rt`, `d_ri`, `read_more`.

Regex literals in the binary (Observed in `strings`): `^([^\d:]+):\s*(.+)` (generic key:value line),
`\(\s*(\d+\.?\d*)\s+(\d+\.?\d*)\s*\)` (peak pair "(mz intensity)"), `(\d+\.?\d*)` (number),
`^\s*(.*?)\s*$` (trim).

Header keys (Observed strings): `NAME:`, `QUANTIFICATION:`, `RESOLUTION:`, `NUM PEAKS:`, plus `NAME_GMB`
and `MST N` (purpose Unknown). The `Initial import from` string suggests a progress or log message.

| Key | Parsed into | Tag |
|---|---|---|
| NAME: | compound name | Observed (string); Inferred (use) |
| QUANTIFICATION: | quant_ions (floats; set) | Observed (string); Inferred (parse) |
| RESOLUTION: | resolution (default "1.0") | Observed (string and default); Inferred |
| NUM PEAKS: | expected peak count | Observed (string); Inferred |
| peak lines `(mz int)` | ions and values | Observed (regex); Inferred (line format) |
| other `key: value` lines | features map (string to string) | Observed (regex); Inferred |

Timing fields: `RT:` is minutes and becomes the spectrum rt in ms (Observed: an MSL written by the old
export and read back by `importMSL` gives the same rt as the source library, `msl_v5.lbr`, 2026-10-20).

A record without `RI:` or `RT:` leaves both at -1 (unknown): the old `importMSL` and `exportMSL` give `RI:-1.0000` and
`RT:-1.0000` for it, while an explicit `RI:0` and `RT:0` stay 0 (Observed, 2026-10-20; the search takes `RI >= 0` as
known, §5.6 in `algorithms.md`).

Rewrite choices (`md_core::library`): `RI:` is the spectrum ri, `RT:` is minutes stored as ms (matches the
old program, checked in `lbr_sample.rs`), `RESOLUTION:` is the spectrum bin width, `NUM PEAKS:` must
equal the number of pairs. Imported compounds get `comb_spec` = spectrum, `ion_freq` 1 per ion and
`spec_count` 1 like the sample library. Export writes the same keys, then features, then 5 pairs per line.

### 6.2 MSL export: `LibrarySearch::exportMSL(path, ExportType)` (Observed; ported in `md_core::library::write_msl_as`)

`DialogExportLibrary::on_okButton_clicked` (0x4cd3b4) maps the combo text to the type: "MSL Nominal" = 1,
"MSL Profile" = 0, "MSL FWHM" = 2, then calls `exportMSL` (0x629f70). The file is opened for writing; failure
throws `IOException` "Error opening library file for writing." and a write error "Error writing data to library
file.". Per compound (in library order; all numbers `std::fixed`):

```
NAME:<name>
RI:<ri, 4 decimals>
RT:<rt in minutes, 4 decimals>
<key>: <value>                     one line per feature, sorted by key; a ';' in the value starts a new "<key>: " line
QUANTIFICATION:<ion, 4 decimals><tab> ...      every ion followed by a tab
RESOLUTION:<bin width, 4 decimals>              not in the FWHM type
NUM PEAKS:<n>
<5 peaks per line, each followed by a space; a line break before peaks 0, 5, 10, ...>
<blank line>
```

and one more blank line at the end of the file. Peak layout by type (`(` m/z, intensity, `)`):

| Type | Peak text | Example |
|---|---|---|
| 1 Nominal | m/z width 4, 0 decimals; intensity width 4 | `(  41   80) ` |
| 0 Profile | m/z width 10, 6 decimals; intensity width 4 | `( 41.000000   80) ` |
| 2 FWHM | m/z width 10, 6 decimals; FWHM width 5, 2 decimals; intensity width 4 with no space before it | `( 41.000000  1.00  80) ` |

Before writing, types 1 and 2 ask the compound for its centroid compound (`LibraryCompound` virtual
`getCentroidLibraryCompound`, 0x61e2ce, which calls `Compound::getCentroidCompound(map&)` 0x611fa8); Nominal then
bins it to 1.0 (`getBinnedCompound`, 0x60dfc8). On an exception the old code prints "ERROR: Error getting
centroid compound for MSL export. Using binned compound." to stderr and writes the binned compound.

- **Centroiding** (only without a FWHM list; with one the spectrum is written as it is): the points are cut into
  groups. A group climbs until two falling steps in a row, then falls
  until a clear rise, and ends at a gap of more than 0.01 in mass. A lone point stays with FWHM 1.00. In a longer
  group the points with at least half of the group maximum count: mass = intensity-weighted mean, intensity = sum,
  FWHM = number of points. A group without weight is dropped (`determineCentroid`).
- **Nominal** (`getBinnedSpectrum`, width 1): the peaks of one whole mass (x.5 rounds up) become one peak with the
  largest intensity (not the sum), and `RESOLUTION` is 1.0000.
- **Intensities** of all three types go through `getNormalizedIntensities`: round(v / base peak x 1000), a result of
  0 becomes 1, so the base peak is 1000.

Checked against the old program (`tools/golden/make-msl-export.sh`, which loads a `.lbr` with
`LibrarySearch::fromDisk` or a `.msl` with `importMSL` and exports all three types): the rewrite writes the same
bytes for all 27 compounds of `lib_aa.lbr` (test `msl_export_matches_the_old_program`, needs the sample files, which
are not committed) and for about 2,250 random compounds with close masses, ties at x.5 and different resolutions
(seeds 1 to 9, checked once by hand, the generator is not committed).

Not ported: line breaks in names and values (the rewrite replaces them by a space so a name cannot start a record;
the old program writes them raw), duplicate masses (the old program throws in `getNormalizedIntensities` and falls
back), a lone point at mass 0 (dropped by the old code) and a base peak of 0. The files are not read back exactly by
`parse_msl` (four decimals) and `parse_msl` does not read the FWHM type. The FWHM type joins FWHM and intensity
without a space when the intensity has four digits (`( 57.000000  1.001000)`): no parser can split that peak, the
old program writes it the same way.

### 6.3 NIST text (`NISTParser`, `SQLLibrary::importNIST`)

Decoded with a black-box driver that calls the old binary's own functions on a synthetic file
(`tools/nist/make-nist-dump.sh`, input `crates/md-core/tests/data/synthetic-nist.sdf`, invented names and
peaks; no NIST data). Reader: `md_core::library::parse_nist` (`.sdf`), tested on the same file
(`tests/nist_synthetic.rs`).

Grammar (Observed, 10 records parsed):
- SDF-like records end with a line `$$$$`.
- A field starts with a line matching `\s*>\s*<(.*)>`; the name is the text between `<` and `>` and may
  contain spaces (`ESTIMATED KOVATS RI`).
- The value is the following lines up to the next field or `$$$$`, each trimmed with `" \n\r\t"`, blank
  lines dropped, joined with `\n`.
- Fields may be missing and come in any order. Lines before the first field (an SDF mol block) and a
  leading UTF-8 byte order mark are ignored by the rewrite (Inferred for the old parser).

Import (`importNIST` 0x6e4024, one SQLite `compound` row per record):

| Field | Into | Absent |
|---|---|---|
| `NAME` | name (a `"` becomes `''` in the SQL text only: SQL escaping) | empty |
| `ESTIMATED KOVATS RI` | RI | -1 |
| `RI ESTIMATION ERROR`, `FORMULA`, `CASNO`, `NISTNO`, `ID`, `MW`, `SYNONYMS` | own columns; rewrite: features under the field name | 0, `''` or -1 |
| `MASS SPECTRAL PEAKS` | one `m/z intensity` pair per line, m/z by strtol, intensity by strtod | no record |

- Only peaks with m/z > 69 and intensity > 0 are kept (Observed: 50 and 60 dropped). A record without a
  kept peak is skipped with `No intensity, skipping...` (Observed).
- RT is -1. Unknown fields: the rewrite keeps them as features (old: not tested).
- Constants `factor_intensity` 0.6 and `factor_mz` 2.0 and the computed `spec_length` are not needed by
  the rewrite (the search uses the `LS_*` settings).

Fragment index (`red_spectrum`, built by `Spectrum::getSpectrumFragments` and `Fragment`, Observed in the
disassembly and Verified against the stored rows):
1. Walk the peaks in mass order. A peak whose next peak is exactly 1 m/z higher and less intense opens or
   continues a run; the first peak where that fails closes the run and belongs to it. A run is kept when
   `100 · sum / total signal >= 1.0`. A peak outside every run is a lone peak; a run below 1 % is dropped
   (its peaks are not lone).
2. Lone peaks directly below a run (M-1, M-2, … without gaps) join it ("hydrogen peaks").
3. The fragment's ion is the most intense mass (the first on ties), its intensity the sum of the group.
   Examples: 71/72 → 71 (112); 70/71/72 with 5/100/12 → 71 (117); 119/120 with 10/100 → none (rising);
   133/134/135 with 100/10/50 → 133 (110), 135 stays alone.

Pre-filter (`SQLLibrary::getPotentialHits` 0x6e2350): the query's fragments are sorted by summed intensity,
strongest first (`compare_intens`), and the first `SQL_MAX_NUMBER_FRAGMENTS` = 5 are matched against the
stored ions with `frag = ion`. A compound is a candidate with `hits > SQL_MINIMAL_NUMBER_EQUAL_FRAGMENTS`
(3), so at least 4 shared fragments (Verified: 4 shared → candidate, 3 shared → not). Both values come
from the settings dialog (§5); the rewrite has them in `SearchSettings::max_fragments` and
`min_equal_fragments`. The tolerance
argument is never read. Query ions below m/z 70 are not dropped: a query whose strongest fragments lie
below 70 (41/42, 43/44, 55/56, 57/58) finds no candidate, because the library side never stores them
(Verified: no candidate and 0 hits in the old binary, the same in the rewrite). The candidates are then
scored like `LibrarySearch` (§5 of algorithms.md):
Verified for an identical spectrum (S = 1, RI 0.969233234, O = 0.989637399) and for a different one
(S = 0.773023519, O = 0.842290311); the rewrite agrees to 1e-8.

Also present: `Name:`, `Formula:`, `Comment:`, `Ions:` strings in the GUI or compound-info text (purpose
Unknown; probably the display/export of a compound, not this reader).

### 6.4 NIST MSP text (public format, not in the old program)

`parse_msl` also reads MSP files as published by NIST and MoNA/MassBank: keys `Name:`, `Synon:`, `Formula:`,
`MW:`, `CAS#:`, `NIST#:`, `Comments:`, `Num Peaks:`; after `Num Peaks:` a line that starts with a digit holds
`m/z intensity` pairs separated by `;` or spaces, each optionally followed by a quoted annotation. Repeated
keys (`Synon:`) are joined with `\n`. Retention-index notations in comments are not read.

### Unknowns (section 6)
- Full MSL line grammar, rt/ri key names, and handling of multiple spectra per block.
- Whether `NUM PEAKS` is trusted or recounted.
- NIST: duplicate m/z in one record, ties in the fragment order. (`LS_IDENTIFICATION_ALGORITHM` 1 and 2
  are unreachable, §5.)
- Encoding of MSL and NIST files.

---

## Unknowns and gaps (summary)

1. .bin/.idx: DataType enum names; compression thresholds (layout decoded, section 2; the TOF polynomial is read since #29, section 4.3).
2. .cmp: versions 1 to 3 (no sample); meaning of `count` and of the model peak's ion code.
3. .calib: encoding only (layout, unit and members verified, section 4.1).
4. .mzcalib: encoding; automatic calibration. FastFlight: time unit, the meaning of `b`, the compressing path.
5. GCMSSettings: fromStream read order and checksum check; persistence of LS_CUTOFF (SQL_* are QSettings).
6. LBR: Peak<int,float> records of version 2/3 compounds with peaks (decoded otherwise, section 1).
7. MSL: most functions undecoded; grammar inferred from strings. NIST: decoded (6.3) except the gaps in its list.
8. Ghidra no-return truncations hide the tail of several readers (fromStream, readIndexFile, readTOFCalibData,
   RICalculator::fromDisk, importMSL). A disassembly-level or Ghidra no-return fix is needed to close these.
