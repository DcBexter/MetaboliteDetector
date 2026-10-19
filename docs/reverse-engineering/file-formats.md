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
The sample file has `version = 2` (produced by an older writer). Parser support: version 2 only.

### 1.1 Top-level fields

| Offset | Type | Meaning | Tag |
|---|---|---|---|
| 0 | u32 | version (writer 5; sample 2) | Observed (toDisk 0x721070); sample value Observed (byte check) |
| 4 | u32 | compound count | Observed (toDisk 0x721070) |
| 8 | record | LibraryCompound records (see 1.2) | Observed |

### 1.2 LibraryCompound record (per compound)

Order as implemented by `parse_lbr.py` (validated on all 27 sample compounds, every check value matches):

| # | Type | Meaning | Tag |
|---|---|---|---|
| 1 | Spectrum | reference spectrum (see 1.4) | Observed (Spectrum::toStream 0x7027a2, verified by parser) |
| 2 | Compound | compound header (see section 3) | Observed (Compound::toStream 0x70fcba, verified by parser) |
| 3 | u32 count + (f32 mz, i32 intensity)[count] | `comb_spec` (combined spectrum) | Observed (LibraryCompound::toStream 0x71f832, element types Inferred from parser) |
| 4 | u32 count + (f32, u32)[count] | `ion_freq` (ion frequency map) | Observed (toStream 0x71f832); element types Inferred |
| 5 | u32 | `spec_count` | Observed (toStream 0x71f832) |
| 6 | u32 count + f64[count] | `ris` (retention index list) | Observed (toStream 0x71f832); element type Inferred |
| 7 | u32 count + 40-byte block records | **unresolved** block list between `ris` and `sources` (see 1.3) | Observed in sample bytes; semantics Unknown |
| 8 | u32 count + string[count] | `sources` | Observed (toStream 0x71f832) |

Not present in the sample data (the parser raises NotImplementedError): `formulas` and `mid_sizes`, which
the decompiled `toStream` lists after `sources`. Their layouts are Unknown.

### 1.3 Block record (40 bytes, sample only)

| Offset in record | Type | Value in sample | Meaning | Tag |
|---|---|---|---|---|
| 0 | f32 | m/z | | Observed (byte check) |
| 4 | f64 | 1.0 | constant | Observed (byte check) |
| 12 | f64 | intensity | equals the comb_spec intensity for the same m/z | Observed (byte check) |
| 20 | f64 | 0.0 | constant | Observed (byte check) |
| 28 | f64 | -1.0 | constant | Observed (byte check) |
| 36 | u32 | 1 | constant | Observed (byte check) |

Discrepancy: the v5 `LibraryCompound::toStream` field list in the decompile does not show a block list
between `ris` and `sources`, but the version 2 sample contains one. The meaning of the block list
(likely a version 2 model or peak list) is an open gap.

### 1.4 Spectrum (used in 1.2 #1)

Observed: `Spectrum<int,float>::toStream` 0x7027a2 (decompile) and byte check in `parse_lbr.py`.

| Offset | Type | Meaning |
|---|---|---|
| 0 | u32 | scan number |
| 4 | f64 | rt field (unit Inferred: ms) |
| 12 | f64 | ri |
| 20 | f32 | bin_width |
| 24 | u32 | n_masses |
| 28 | f32[n_masses] | masses |
| ... | u32 | n_values |
| ... | i32[n_values] | values (intensities) |
| ... | u32 | check = (scan + n_masses) mod 2^32 (Observed in sample) |

### 1.5 Validated parser

`docs/reverse-engineering/tools/parse_lbr.py`. Run with `python3 -I parse_lbr.py < file.lbr` (or a decompressed `.raw`).
Output on the sample: version 2, 27 compounds, all bytes consumed, first compound Hexatriacontane (ri 3599.987).
Accepts gzip input automatically. Only version 2 is decoded; version 5 raises NotImplementedError.
The Rust reader `md_core::library::parse_lbr` returns the same values for every field of all 27 sample
compounds (`tools/dump_lbr.py` against `md-core/tests/lbr_sample.rs`, 2026-10-19). It rejects version 5
and compounds with peak records as not supported.

### Unknowns (section 1)
- Meaning of the 40-byte block list (1.3) and its relation to the v5 layout.
- `formulas` and `mid_sizes` layouts (absent in the sample).
- Version 5 peak/model layout (`Peak<double,short>` and `PeakEssential<short>`, see section 3.3).
- Element types of comb_spec/ion_freq/ris (inferred from the parser, not from DWARF).
- String encoding (latin-1 assumed by the parser; the real encoding is unverified).

---

## 2. Chromatogram data (.bin + .idx)

Symbols: `GCMSDiskScanWriter<int,float>` (writer), `GCMSDiskScan<int,float>` (reader, wraps a
`DataProvider<int,float>`), `DataProvider<int,float>::readIndexFile` / `writeIndexFile`.

File naming
- Index file = data file name + `".idx"` (Observed: `DataProvider::readIndexFile` 0x75bf30 appends ".idx"
  to `f_name` and opens it `in|binary`; `GCMSDiskScanWriter::writeIndexFile` 0x6e7d00 appends ".idx" to
  `cache_file` and opens `out|trunc|binary`).
- Data file = the name passed to `GCMSDiskScanWriter(char const*, DataType)` 0x6e610c, stored as the
  cache file name. The extension `.bin` is Inferred from the brief, not verified in code.

### 2.1 Scan record (data file)

Written by `GCMSDiskScanWriter::addScan(uint scan, double time, vector<float>& masses,
vector<int>& values, streampos pos, double tic)` 0x6e6888 (Observed: the decompiled write sequence below).
The vector overload 0x6e7072 calls it with `pos = 0` and `tic = -1.0` (Observed), so that path writes the
masses. The other overloads (0x6e6fd6, 0x6e65c8) were not decompiled.

Write order (Observed; sizes from the `write` argument types and the `size_all` arithmetic):

| Offset in record | Type | Meaning | Tag |
|---|---|---|---|
| 0 | i32 | `size_all`: total record length in bytes including this field. Masses written: 36 + 8n. Masses shared: 36 + 4n | Observed (0x6e6888) |
| 4 | i64 | `masses_ptr`: absolute file offset of the mass array (record start + 36 + 4n when masses written, else the given pos) | Observed (0x6e6888) |
| 12 | f64 | `tic`: total ion current; computed by `calculateTIC` when passed -1.0 | Observed (0x6e6888) |
| 20 | u32 | scan number | Observed (0x6e6888) |
| 24 | f64 | retention time (`time`) | Observed (0x6e6888) |
| 32 | u32 | n = number of values (ions) in this scan | Observed (0x6e6888) |
| 36 | i32[n] | intensities (`s_values`) | Observed (0x6e6888) |
| 36 + 4n | f32[n] | masses (`s_masses`), present only when masses are written | Observed (0x6e6888) |

Header of the data file: none observed. The first record starts at offset 0 (Inferred: addScan begins
with `seekp(..., _S_beg)` whose offset operand Ghidra did not recover, so the start position is Unknown).

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
- `CENTROID` type (Observed, 0x6e6888): stores a zeroed BaselineStatistics and writes the data unchanged.
- Other profile types: data written unchanged (Observed, 0x6e6888: no branch).

### 2.3 Index file (.idx)

Writer: `GCMSDiskScanWriter::writeIndexFile` 0x6e7d00. Reader: `DataProvider::readIndexFile` 0x75bf30.
Both decompiles are **truncated** right after opening the file (`c_str` no-return artifact). The field list
below is **Inferred** from the local variable names in both functions, not from the byte-level writes.

| Item (inferred order) | Inferred type / content | Evidence (local names) |
|---|---|---|
| magic number | i64 (value not recovered) | `magic_number` (longlong) in both functions |
| type size markers | u32 values for sizeof(int), sizeof(double), sizeof(long long) | `size_int`, `size_double`, `size_long_long` |
| TOF calibration polynomial | u32 `poly_size` + f64[poly_size] | `poly_size`; `setTOFCalibrationPolynomial` 0x6e71b6 |
| all ions (masses) | u32 `size_all_ions` + f32[] | `size_all_ions`, `ions` set |
| scan count | u32 | `GCMSDiskScan::initialization` 0x6e5ba6 uses `getScanCount` |
| per-scan entries | scan (u32), rt (f64), tic (f64), record pos (i64), mass pointer (i64), num_val (u32), i_compressed (int) | `scan`, `rt`, `tic`, `pos`, `mass_ptr`, `num_val`, `i_compressed` |
| baseline statistics per compressed scan | mean, sd, msd, nf (f64 each) | `mean`, `sd`, `msd`, `nf` |
| features | u32 count + pairs of length-prefixed strings | `f_count`, `s_key`, `s_value`, `key`, `value` |
| RI calibration | embedded `RICalculator` (see section 4.1) | local `RICalculator ric` in `writeIndexFile` |

### 2.4 Reading

`GCMSDiskScan::getSpectrum(unsigned int)` 0x6e2638 (truncated): calls `DataProvider::moveToScan`,
`getTime`, `getIonCount`, `getSignals(0, n-1)`, `getMasses(0, n-1)`. The provider implementation that
seeks into the data file using the index tables is not decoded (gap).

### Unknowns (section 2)
- Byte-level index layout and magic value (writer and reader are truncated; the magic is not recovered).
- `DataType` enum values (PROFILE, CENTROID, COMPRESSED_PROFILE, ...).
- `setCompressionParameter` default thresholds and the BaselineStatistics field order.
- Meaning of the seekp offset in `addScan` (start position of the record).
- Whether the data file has a header (none observed in the writer).
- The DataProvider reader (seek logic) and `GCMSDiskScanWriter::close`.
- File extension `.bin` (not verified in code).

---

## 3. Detected compounds (.cmp)

Symbols: `Compound<int,float>` (header), `Peak<int,float>`, `PeakEssential<float>`, `Peak<double,short>`,
`PeakEssential<short>`.

### 3.1 File container

- `Compound<int,float>::toDisk(vector<Compound*>, file)` 0x711332 (Observed): opens an `ogzstream`
  (gzip), calls `compoundsToStream`, closes. So `.cmp` is gzip-compressed.
- `compoundsToStream(vector<const Compound*>, ostream)` 0x7111d4 (Observed): `u32 version = 4`, `u32 count`,
  then `Compound::toStream` for each entry. The version is 4 for `.cmp` (the library writer uses 5).

| Offset | Type | Meaning | Tag |
|---|---|---|---|
| 0 | u32 | version = 4 | Observed (0x7111d4) |
| 4 | u32 | compound count | Observed (0x7111d4) |
| 8 | Compound record x count | | Observed (0x7111d4) |

Checksum: there is no file-level checksum. Each compound record ends with a u32 check value (3.2).

### 3.2 Compound record (Compound<int,float>::toStream 0x70fcba)

Observed in the earlier decompile of `toStream` (0x70fcba) and cross-checked against all 27 sample
compounds by `parse_lbr.py` (`read_compound`). The same record appears inside LibraryCompound (section 1.2 #2).

| # | Type | Meaning | Tag |
|---|---|---|---|
| 1 | u32 | n_peaks (count of PeakEssential<float> records that follow; 0 in the sample) | Observed (sample), Inferred (layout of n>0) |
| 2 | string | name | Observed (sample) |
| 3 | i32 count + f64[count] | fwhm list; count -1 = null | Observed (sample) |
| 4 | i32 count + f64[count] | integrals list; count -1 = null | Observed (sample) |
| 5 | u32 count + f32[count] | quant_ions (`std::set<float>`, sorted) | Observed (sample) |
| 6 | u32 count + f32[count] | auto_quant_ions (`std::set<float>`) | Observed (sample) |
| 7 | u32 count + (string, string)[count] | features (`std::map<string,string>`) | Observed (sample) |
| 8 | f64 | drt (`Compound::drt`, offset 376 in DWARF) | Observed (sample), name from DWARF |
| 9 | i64 | time_stamp (`time_t`, DWARF offset 360) | Observed (sample), name from DWARF |
| 10 | i64 | id (DWARF offset 368) | Observed (sample), name from DWARF |
| 11 | u32 | check = (len(quant_ions) + len(features) + len(name bytes) + n_peaks) mod 2^32 | Observed (sample, all 27 match) |

Note: `Compound` members (DWARF): name@120 std::string, features@312 map, time_stamp@360, id@368, drt@376.

### 3.3 Peak records (only when n_peaks > 0)

- `Peak<double,short>::toStream` 0x6f66e2 (Observed, earlier pass): `PeakEssential<short>::toStream`, then
  `i32 nf`, `u32 n`, `f64[n]`, then `u32 check = n + getFirstScan + getLastScan + getMaxScan`.
- `PeakEssential<short>::toStream` 0x6facfe (Observed, earlier pass): fields in order `ion`, `sharpness`,
  `integral`, `quality`, `sn`, `count` (uint), `begin`, `end`, `model_parameter`, `max_scan`,
  `interpolated_max_scan`. Types were not re-checked in this pass.
- `PeakEssential<float>::toStream` (nm 0x5f994c, Ghidra 0x6f994c) and `Peak<int,float>::toStream` (nm 0x5f3036,
  Ghidra 0x6f3036): **not decoded**. The `.cmp` peak layout for float peaks is Unknown.

### 3.4 Reader

`Compound<int,float>::fromStream(istream&, size_t* version)` 0x7103cc: Ghidra decompile truncated at
`operator new(0x180)` (Observed: object size 0x180 bytes). `Compound::fromDisk` 0x711896 and
`compoundsFromStream` 0x711596 are not decoded. The reader order is Inferred to mirror 3.1 and 3.2.

### Unknowns (section 3)
- Layout of PeakEssential<float> and Peak<int,float> (needed for n_peaks > 0).
- Whether the reader checks the version or accepts only version 4.
- Fields between fwhm and features in the reader if they differ from the writer (not recovered).

---

## 4. Calibration files

### 4.1 RI calibration (.calib): `gcms::RICalculator`

Members: `rt2ri` (`std::map<double,double>`), `drt` (f64), and a `std::string` (name not recovered;
DWARF lookup not completed).

Binary form, `RICalculator::toStream` 0x7d2fd2 (Observed):

| Offset | Type | Meaning |
|---|---|---|
| 0 | u32 | entry count (`rt2ri.size()`) |
| 4 | f64 rt, f64 ri | repeated count times, in ascending rt order (std::map iteration) |
| 4 + 16n | f64 | drt |
| ... | u32 | string length |
| ... | bytes | string data (raw, no terminator) |

Text form (the `.calib` file), `RICalculator::toDisk(file)` 0x7d28f0 (Observed for format, Inferred for loop):
- Opens the file with `ofstream(trunc)`, writes each entry as `<rt>` TAB `<ri>` followed by `endl`.
- Values are written with the default ostream precision (6 significant digits, Observed: no `setprecision`).
  The text form is therefore lossy compared with the binary form.
- The decompile shows only one iteration because Ghidra cut the loop after the no-return `endl` call (Inferred).
- Throws if the file cannot be opened (Observed: `__cxa_allocate_exception` on failure).

Reader, `RICalculator::fromDisk(file)` 0x7d22ac: clears `rt2ri`, opens the file, reads the first line
with `getline`, and splits it with `boost::regex("\\s+")` into rt and ri (Observed). The rest of the loop
is truncated (Inferred: one line per entry, whitespace-separated, blank lines handled by the regex).
Note: `toDisk` writes TAB, `fromDisk` splits on any whitespace, so both work.

Field table for the text file:

| Column | Type (text) | Meaning |
|---|---|---|
| 1 | decimal (6 significant digits) | rt (unit Inferred: min or ms, see unknowns) |
| 2 | decimal (6 significant digits) | ri |

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

- Precision of the numbers: `QTextStream::setRealNumberPrecision` is called, but its argument was not recovered
  (Unknown; Qt default is 6 significant digits).
- Encoding: `QTextStream` default (Unknown; ASCII content only in the strings).

Reader: `MainWindow::readTOFCalibData(QString file, QWidget*)` 0x4c8aba. Decompile truncated at
`QString::toStdString` / `std::string::c_str`. Locals (`order`, `tmp`, `i`, `rc_par`, `char text[32]`)
suggest a line-based `ifstream` reader that parses the keys above (Inferred).

`TOFCalibrationWindow::saveSettings()` 0x6c077e is **not** the .mzcalib writer. It stores only the window
size in `QSettings("TU-Braunschweig", "TOFCalib")`, group "GUI", key "size" (Observed).

### Unknowns (section 4)
- RI unit of rt in `.calib` (min or ms; Inferred from the GUI label "Retention time of displayed spectrum (ms)"
  and the RI table labels, not verified).
- Name and meaning of the `std::string` member of RICalculator (DWARF lookup not completed).
- Exact precision in `.mzcalib` (setRealNumberPrecision argument not recovered).
- Encoding of both text files.
- The complete reader for `.mzcalib` and `.calib` (both truncated).

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
- Persistence of LS_CUTOFF, AN_MODELED_VALUES_FOR_INTEGRATION and the SQL_* fields (not in this stream).
- Exact types of LS_RESOLUTION (f32 inferred) and LS_IDENTIFICATION_ALGORITHM (enum width 4 inferred).
- Whether the AN_MAX_BASE_DEV static write in GCMSSettingsObject is intentional.

---

## 6. MSL (AMDIS) and NIST MSP import/export

Status: partially decoded. The import function is decompiled but truncated; the NIST parser functions are
named but not decoded. Most of this section is based on string literals and local names (Inferred).

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

Timing fields (`rt`, `ri`): the importer stores them as `d_rt` and `d_ri`. Which key holds them is not
recovered (Gap).

Rewrite choices (`md_core::library`, not verified against the old app): `RI:` is the spectrum ri, `RT:` is
minutes (AMDIS convention) and is stored as ms, `RESOLUTION:` is the spectrum bin width, `NUM PEAKS:` must
equal the number of pairs. Imported compounds get `comb_spec` = spectrum, `ion_freq` 1 per ion and
`spec_count` 1 like the sample library. Export writes the same keys, then features, then 5 pairs per line.

### 6.2 MSL export: `DialogExportLibrary` (not decompiled)

Strings: `Error getting centroid compound for MSL export. Using binned compound.` and the same header keys as
6.1 (`NAME:`, `QUANTIFICATION:`, `RESOLUTION:`, `NUM PEAKS:`). The writer order is Inferred to match 6.1.
Gap: the export function was not decoded.

### 6.3 NIST MSP parser (`NISTParser`)

Named functions (addresses from the earlier nm pass, Ghidra = nm + 0x100000): constructor 0x7ed026,
nextRecord 0x7ed2be, nextField 0x7ed448, getFieldName 0x7ed88a, getFieldValue 0x7ed94c, trim 0x7edabe.
**None of these were decompiled in this pass.** Gap.

Field keys present as strings (Observed; their use by the `ImportNist` SQL path is Inferred from string adjacency):
`NAME`, `SYNONYMS`, `FORMULA`, `CASNO`, `NISTNO`, `ESTIMATED KOVATS RI`, `RI ESTIMATION ERROR`,
`MASS SPECTRAL PEAKS`. The SQL insert writes
`compound (time_stamp, id, name, ri, ri_err, synonyms, formula, mw, cas, nist_no, nist_id, factor_intensity, ...)`.

Also present: `Name:`, `Formula:`, `Comment:`, `Ions:` strings in the GUI or compound-info text (purpose
Unknown; probably the display/export of a compound, not the MSP reader).

### Unknowns (section 6)
- Full MSL line grammar, rt/ri key names, and handling of multiple spectra per block.
- Whether `NUM PEAKS` is trusted or recounted.
- NIST parser state machine and field order (all functions undecoded).
- Encoding of MSL and NIST files.
- Export order and number formatting for MSL.

---

## Unknowns and gaps (summary)

1. .bin/.idx: index byte layout and magic value; DataType enum; compression thresholds; reader seek logic.
2. .cmp: PeakEssential<float> and Peak<int,float> layouts; reader details.
3. .calib: rt unit; RICalculator string member; loop structure of toDisk (Inferred).
4. .mzcalib: reader; precision; encoding.
5. GCMSSettings: fromStream read order and checksum check; persistence of LS_CUTOFF and SQL_* fields.
6. LBR: block list meaning (1.3); formulas and mid_sizes; version 5 peak/model layout.
7. MSL and NIST: most functions undecoded; grammar inferred from strings.
8. Ghidra no-return truncations hide the tail of several readers (fromStream, readIndexFile, readTOFCalibData,
   RICalculator::fromDisk, importMSL). A disassembly-level or Ghidra no-return fix is needed to close these.
