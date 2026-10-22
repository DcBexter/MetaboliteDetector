//! The old `GCMSSettings` (file-formats.md §5): all values in one struct, the old file format, and the
//! conversion to the settings of each analysis step. Files are untrusted: `from_bytes` checks size,
//! checksum and every value before anything is used.

use std::io::Read;

use flate2::read::GzDecoder;

use crate::batch::{QuantSettings, ScanSettings};
use crate::deconv::DeconvSettings;
use crate::peaks::{Detector, PeakSettings};
use crate::ri::RiCalibration;
use crate::search::SearchSettings;

/// The settings of every analysis step, derived once from [`GcmsSettings`] by [`GcmsSettings::analysis`].
/// The batch scan, the analyze run, the md-cli and the app's batch window read this one struct.
#[derive(Debug, Clone, PartialEq)]
pub struct AnalysisSettings {
    pub peaks: PeakSettings,
    pub deconv: DeconvSettings,
    pub quant: QuantSettings,
    pub search: SearchSettings,
    /// `SF_FWHM`: share of a peak's maximum a point must reach to count in its centroid. Only
    /// [`IonTraces::new`] reads it; a change re-runs the analysis (the app keys its analysis on it).
    pub fwhm: f64,
}

impl Default for AnalysisSettings {
    /// The step defaults of the old static values, without the settings file's automatic width.
    fn default() -> Self {
        Self {
            peaks: PeakSettings::default(),
            deconv: DeconvSettings::default(),
            quant: QuantSettings::default(),
            search: SearchSettings::default(),
            fwhm: 0.5,
        }
    }
}

/// Old `GCMSSettings`. Defaults are the old static values (file-formats.md §5.1).
/// Fields marked "stored only" are kept for the file but no step of the rewrite reads them yet.
#[derive(Debug, Clone, PartialEq)]
pub struct GcmsSettings {
    // Spectrum filter. In 3.3 the first three are only stored (no analysis reads them), `SF_FWHM` steers
    // every centroiding and `SF_CENTROID_*` the centroid export and FastFlight import (algorithms.md §4b.4).
    pub sf_peak_threshold_begin: f64,
    pub sf_peak_threshold_end: f64,
    pub sf_max_base_dev: f64,
    pub sf_fwhm: f64,
    pub sf_centroid_threshold_begin: f64,
    pub sf_centroid_threshold_end: f64,
    pub sf_centroid_max_base_dev: f64,
    // Analysis.
    pub an_peak_threshold_begin: f64,
    pub an_peak_threshold_end: f64,
    /// Stored only: the old peak detector does not read it (algorithms.md §1).
    pub an_max_base_dev: f64,
    pub an_min_peak_height: f64,
    /// `None`: automatic, from each run's scan rate (`deconv::auto_width`). The old file holds a number,
    /// so a loaded file sets it.
    pub an_deconvolution_width: Option<f64>,
    pub an_number_bins: i32,
    pub an_min_peak_number: i32,
    // Library search.
    pub ls_ri_diff: f64,
    pub ls_pure_factor: f64,
    pub ls_impure_factor: f64,
    pub ls_scale_factor_ion: f64,
    pub ls_scale_factor_value: f64,
    pub ls_scale_factor_no_match: f64,
    /// Stored only.
    pub ls_resolution: f32,
    pub ls_use_ri: bool,
    /// Stored only.
    pub ls_use_scaled_lib: bool,
    pub ls_use_spec: bool,
    // Quantification (ion number, distance and cutoff are stored only; excluded ions steer the model peak).
    pub quant_ion_number: i32,
    pub quant_min_dist_ions: f64,
    pub quant_quality_cutoff: f64,
    pub quant_excl_ions: Vec<(f32, f32)>,
    // Not in the old file: they keep their value when a file is loaded.
    /// `AN_MIN_HEIGHT_BP`.
    pub an_min_height_bp: f64,
    /// `AN_MODELED_VALUES_FOR_INTEGRATION`.
    pub an_modeled_values_for_integration: bool,
    /// `LS_MASS_FILTER`.
    pub ls_mass_filter: Vec<(f32, f32)>,
    /// `LS_CUTOFF`: a hit below this overall score counts as no match (result table, Find substance, batch).
    pub ls_cutoff: f64,
    /// `SQL_MAX_NUMBER_FRAGMENTS` and `SQL_MINIMAL_NUMBER_EQUAL_FRAGMENTS`: NIST pre-filter.
    pub sql_max_number_fragments: i32,
    pub sql_minimal_number_equal_fragments: i32,
    /// Peak detection method (old Calibration Setter window, not saved by the old app).
    pub peak_detector: Detector,
    /// QSettings `[Algorithm] baseline_correction` of the main window (algorithms.md §4b.4):
    /// [`crate::Chromatogram::set_baseline_correction`] on every opened run. Batch runs stay uncorrected.
    pub baseline_correction: bool,
}

impl Default for GcmsSettings {
    fn default() -> Self {
        Self {
            sf_peak_threshold_begin: 10.0,
            sf_peak_threshold_end: -10.0,
            sf_max_base_dev: 20.0,
            sf_fwhm: 0.5,
            sf_centroid_threshold_begin: 10.0,
            sf_centroid_threshold_end: -5.0,
            sf_centroid_max_base_dev: 30.0,
            an_peak_threshold_begin: 10.0,
            an_peak_threshold_end: -5.0,
            an_max_base_dev: 20.0,
            an_min_peak_height: 20.0,
            an_deconvolution_width: None,
            an_number_bins: 10,
            an_min_peak_number: 10,
            ls_ri_diff: 20.0,
            ls_pure_factor: 0.5,
            ls_impure_factor: 0.5,
            ls_scale_factor_ion: 2.0,
            ls_scale_factor_value: 0.6,
            ls_scale_factor_no_match: 1.0,
            ls_resolution: 0.3,
            ls_use_ri: true,
            ls_use_scaled_lib: true,
            ls_use_spec: true,
            quant_ion_number: 3,
            quant_min_dist_ions: 0.5,
            quant_quality_cutoff: 1.0,
            quant_excl_ions: Vec::new(),
            an_min_height_bp: 0.0,
            an_modeled_values_for_integration: true,
            ls_mass_filter: Vec::new(),
            ls_cutoff: 0.8,
            sql_max_number_fragments: 5,
            sql_minimal_number_equal_fragments: 3,
            peak_detector: Detector::Default,
            baseline_correction: false,
        }
    }
}

/// Old `toStream` size without excluded ions, and per excluded ion (§5.2).
const FIXED_LEN: usize = 191;
const ION_LEN: usize = 8;
/// Position of the excluded-ion count; everything before it is fixed-size.
const COUNT_AT: usize = 179;
/// Upper limit on a decompressed file; a real one is under 1 MB even with thousands of ions.
pub const MAX_LEN: u64 = 1 << 20;

/// Reads `f64`/`f32`/`i32`/`u8` in the old little-endian layout.
struct Cursor<'a>(&'a [u8]);

impl Cursor<'_> {
    fn take<const N: usize>(&mut self) -> [u8; N] {
        let (head, rest) = self.0.split_at(N);
        self.0 = rest;
        head.try_into().expect("split_at gave N bytes")
    }
    fn f64(&mut self) -> f64 {
        f64::from_le_bytes(self.take())
    }
    fn f32(&mut self) -> f32 {
        f32::from_le_bytes(self.take())
    }
    fn i32(&mut self) -> i32 {
        i32::from_le_bytes(self.take())
    }
    fn bool(&mut self) -> bool {
        self.take::<1>()[0] != 0
    }
}

/// Bounds of a number without a range of its own.
pub const ANY: (f64, f64) = (-1e9, 1e9);

/// A field of [`GcmsSettings`] as read from and written to the app format, and (for the fixed-size
/// ones) the old file.
pub enum Field<'a> {
    F64(&'a mut f64),
    /// A number, or `auto` for `None`.
    Auto(&'a mut Option<f64>),
    F32(&'a mut f32),
    I32(&'a mut i32),
    Bool(&'a mut bool),
    Ranges(&'a mut Vec<(f32, f32)>),
    Detector(&'a mut Detector),
}

/// The settings keys: the old file and `QSettings` names. `KEYS` and the app's settings window use these,
/// so a typo is a compile error.
pub mod keys {
    pub const SF_PEAK_THRESHOLD_BEGIN: &str = "SF_PEAK_THRESHOLD_BEGIN";
    pub const SF_PEAK_THRESHOLD_END: &str = "SF_PEAK_THRESHOLD_END";
    pub const SF_MAX_BASE_DEV: &str = "SF_MAX_BASE_DEV";
    pub const SF_FWHM: &str = "SF_FWHM";
    pub const SF_CENTROID_THRESHOLD_BEGIN: &str = "SF_CENTROID_THRESHOLD_BEGIN";
    pub const SF_CENTROID_THRESHOLD_END: &str = "SF_CENTROID_THRESHOLD_END";
    pub const SF_CENTROID_MAX_BASE_DEV: &str = "SF_CENTROID_MAX_BASE_DEV";
    pub const AN_PEAK_THRESHOLD_BEGIN: &str = "AN_PEAK_THRESHOLD_BEGIN";
    pub const AN_PEAK_THRESHOLD_END: &str = "AN_PEAK_THRESHOLD_END";
    pub const AN_MAX_BASE_DEV: &str = "AN_MAX_BASE_DEV";
    pub const AN_MIN_PEAK_HEIGHT: &str = "AN_MIN_PEAK_HEIGHT";
    pub const AN_DECONVOLUTION_WIDTH: &str = "AN_DECONVOLUTION_WIDTH";
    pub const AN_NUMBER_BINS: &str = "AN_NUMBER_BINS";
    pub const AN_MIN_PEAK_NUMBER: &str = "AN_MIN_PEAK_NUMBER";
    pub const AN_MIN_HEIGHT_BP: &str = "AN_MIN_HEIGHT_BP";
    pub const AN_MODELED_VALUES_FOR_INTEGRATION: &str = "AN_MODELED_VALUES_FOR_INTEGRATION";
    pub const BASELINE_CORRECTION: &str = "baseline_correction";
    /// Old names of the settings file, kept for the app's alias table (`persist.rs`).
    pub const LEGACY_AN_BP_HEIGHT: &str = "AN_BP_HEIGHT";
    pub const LEGACY_MASS_FILTER: &str = "MASS_FILTER";
    pub const LS_RI_DIFF: &str = "LS_RI_DIFF";
    pub const LS_PURE_FACTOR: &str = "LS_PURE_FACTOR";
    pub const LS_IMPURE_FACTOR: &str = "LS_IMPURE_FACTOR";
    pub const LS_SCALE_FACTOR_ION: &str = "LS_SCALE_FACTOR_ION";
    pub const LS_SCALE_FACTOR_VALUE: &str = "LS_SCALE_FACTOR_VALUE";
    pub const LS_SCALE_FACTOR_NO_MATCH: &str = "LS_SCALE_FACTOR_NO_MATCH";
    pub const LS_RESOLUTION: &str = "LS_RESOLUTION";
    pub const LS_USE_RI: &str = "LS_USE_RI";
    pub const LS_USE_SCALED_LIB: &str = "LS_USE_SCALED_LIB";
    pub const LS_USE_SPEC: &str = "LS_USE_SPEC";
    pub const LS_MASS_FILTER: &str = "LS_MASS_FILTER";
    pub const LS_CUTOFF: &str = "LS_CUTOFF";
    pub const QUANT_ION_NUMBER: &str = "QUANT_ION_NUMBER";
    pub const QUANT_MIN_DIST_IONS: &str = "QUANT_MIN_DIST_IONS";
    pub const QUANT_QUALITY_CUTOFF: &str = "QUANT_QUALITY_CUTOFF";
    pub const QUANT_EXCL_IONS: &str = "QUANT_EXCL_IONS";
    pub const SQL_MAX_NUMBER_FRAGMENTS: &str = "SQL_MAX_NUMBER_FRAGMENTS";
    pub const SQL_MINIMAL_NUMBER_EQUAL_FRAGMENTS: &str = "SQL_MINIMAL_NUMBER_EQUAL_FRAGMENTS";
    pub const PEAK_DETECTOR: &str = "PEAK_DETECTOR";
}

/// One setting: its old key, where it lives in the struct and what values it accepts.
pub struct Key {
    pub name: &'static str,
    /// Display text of the settings window: (English, German).
    pub label: (&'static str, &'static str),
    pub field: fn(&mut GcmsSettings) -> Field<'_>,
    /// Inclusive bounds (m/z for the range lists, the width for the matched detector); `None`: any
    /// number within ±1e9. `validate`, the app's loader and the settings window all read this.
    pub range: Option<(f64, f64)>,
    /// A fixed-size slot of the old file (§5.2); the excluded ions follow them, the rest is not in it.
    pub fixed: bool,
}

const fn key(
    name: &'static str,
    label: (&'static str, &'static str),
    field: fn(&mut GcmsSettings) -> Field<'_>,
    range: Option<(f64, f64)>,
    fixed: bool,
) -> Key {
    Key {
        name,
        label,
        field,
        range,
        fixed,
    }
}

/// Every setting under its old key, in the order of the old file. Parse, both file formats, `validate`
/// and the settings window all go through this one table; the defaults are `GcmsSettings::default()`.
/// The old `QSettings` name `baseline_correction` stays as it is so an old settings file imports it.
// Rustfmt would put each of the 36 rows on five lines.
#[rustfmt::skip]
pub const KEYS: [Key; 36] = [
    key(keys::SF_PEAK_THRESHOLD_BEGIN, ("Spectrum filter: peak start", "Spektrumfilter: Peakbeginn"), |s| Field::F64(&mut s.sf_peak_threshold_begin), None, true),
    key(keys::SF_PEAK_THRESHOLD_END, ("Spectrum filter: peak end", "Spektrumfilter: Peakende"), |s| Field::F64(&mut s.sf_peak_threshold_end), None, true),
    key(keys::SF_MAX_BASE_DEV, ("Spectrum filter: baseline deviation", "Spektrumfilter: Basislinienabweichung"), |s| Field::F64(&mut s.sf_max_base_dev), None, true),
    key(keys::SF_FWHM, ("Centroid: FWHM share", "Zentroid: FWHM-Anteil"), |s| Field::F64(&mut s.sf_fwhm), Some((0.0, 1.0)), true),
    key(keys::SF_CENTROID_THRESHOLD_BEGIN, ("Centroid: peak start", "Zentroid: Peakbeginn"), |s| Field::F64(&mut s.sf_centroid_threshold_begin), None, true),
    key(keys::SF_CENTROID_THRESHOLD_END, ("Centroid: peak end", "Zentroid: Peakende"), |s| Field::F64(&mut s.sf_centroid_threshold_end), None, true),
    key(keys::SF_CENTROID_MAX_BASE_DEV, ("Centroid: baseline deviation", "Zentroid: Basislinienabweichung"), |s| Field::F64(&mut s.sf_centroid_max_base_dev), None, true),
    key(keys::AN_PEAK_THRESHOLD_BEGIN, ("Peak start slope", "Steigung am Peakbeginn"), |s| Field::F64(&mut s.an_peak_threshold_begin), None, true),
    key(keys::AN_PEAK_THRESHOLD_END, ("Peak end slope", "Steigung am Peakende"), |s| Field::F64(&mut s.an_peak_threshold_end), None, true),
    key(keys::AN_MAX_BASE_DEV, ("Peak baseline deviation", "Peak-Basislinienabweichung"), |s| Field::F64(&mut s.an_max_base_dev), None, true),
    key(keys::AN_MIN_PEAK_HEIGHT, ("Minimum S/N", "Minimum-S/N"), |s| Field::F64(&mut s.an_min_peak_height), Some((0.0, 1000000000.0)), true),
    key(keys::AN_DECONVOLUTION_WIDTH, ("Deconvolution width", "Dekonvolutionsbreite"), |s| Field::Auto(&mut s.an_deconvolution_width), Some((0.1, 100.0)), true),
    key(keys::AN_NUMBER_BINS, ("Apex steps per scan", "Apex-Schritte pro Scan"), |s| Field::I32(&mut s.an_number_bins), Some((1.0, 100.0)), true),
    key(keys::AN_MIN_PEAK_NUMBER, ("Minimum ions per component", "Mindestzahl Ionen pro Komponente"), |s| Field::I32(&mut s.an_min_peak_number), Some((1.0, 10000.0)), true),
    key(keys::AN_MIN_HEIGHT_BP, ("Minimum share of base peak", "Mindestanteil am Basispeak"), |s| Field::F64(&mut s.an_min_height_bp), Some((0.0, 1.0)), false),
    key(keys::AN_MODELED_VALUES_FOR_INTEGRATION, ("Integrate fitted model", "Angepasstes Modell integrieren"), |s| Field::Bool(&mut s.an_modeled_values_for_integration), None, false),
    key(keys::BASELINE_CORRECTION, ("Baseline correction", "Basislinienkorrektur"), |s| Field::Bool(&mut s.baseline_correction), None, false),
    key(keys::LS_RI_DIFF, ("RI tolerance (sigma)", "RI-Toleranz (Sigma)"), |s| Field::F64(&mut s.ls_ri_diff), Some((1e-06, 1000000000.0)), true),
    key(keys::LS_PURE_FACTOR, ("Pure spectrum weight", "Gewicht reines Spektrum"), |s| Field::F64(&mut s.ls_pure_factor), Some((0.0, 1000000000.0)), true),
    key(keys::LS_IMPURE_FACTOR, ("Impure spectrum weight", "Gewicht unreines Spektrum"), |s| Field::F64(&mut s.ls_impure_factor), Some((0.0, 1000000000.0)), true),
    key(keys::LS_SCALE_FACTOR_ION, ("Mass exponent", "Massenexponent"), |s| Field::F64(&mut s.ls_scale_factor_ion), None, true),
    key(keys::LS_SCALE_FACTOR_VALUE, ("Intensity exponent", "Intensitätsexponent"), |s| Field::F64(&mut s.ls_scale_factor_value), None, true),
    key(keys::LS_SCALE_FACTOR_NO_MATCH, ("Unmatched peak factor", "Faktor für Peaks ohne Partner"), |s| Field::F64(&mut s.ls_scale_factor_no_match), None, true),
    key(keys::LS_RESOLUTION, ("Search resolution", "Suchauflösung"), |s| Field::F32(&mut s.ls_resolution), None, true),
    key(keys::LS_USE_RI, ("Score", "Score"), |s| Field::Bool(&mut s.ls_use_ri), None, true),
    key(keys::LS_USE_SCALED_LIB, ("Use scaled library", "Skalierte Bibliothek verwenden"), |s| Field::Bool(&mut s.ls_use_scaled_lib), None, true),
    key(keys::LS_USE_SPEC, ("Score", "Score"), |s| Field::Bool(&mut s.ls_use_spec), None, true),
    key(keys::LS_MASS_FILTER, ("Ignored masses", "Ignorierte Massen"), |s| Field::Ranges(&mut s.ls_mass_filter), Some((0.0, 100000.0)), false),
    key(keys::LS_CUTOFF, ("Cutoff score", "Mindest-Score"), |s| Field::F64(&mut s.ls_cutoff), Some((0.0, 1.0)), false),
    key(keys::QUANT_ION_NUMBER, ("Quant ions", "Quant-Ionen"), |s| Field::I32(&mut s.quant_ion_number), Some((0.0, 10000.0)), true),
    key(keys::QUANT_MIN_DIST_IONS, ("Quant ion minimum distance", "Mindestabstand der Quant-Ionen"), |s| Field::F64(&mut s.quant_min_dist_ions), None, true),
    key(keys::QUANT_QUALITY_CUTOFF, ("Quant quality cutoff", "Quant-Qualitätsschwelle"), |s| Field::F64(&mut s.quant_quality_cutoff), None, true),
    key(keys::QUANT_EXCL_IONS, ("Ions left out of the model peak", "Vom Modellpeak ausgenommene Ionen"), |s| Field::Ranges(&mut s.quant_excl_ions), Some((0.0, 100000.0)), false),
    key(keys::SQL_MAX_NUMBER_FRAGMENTS, ("Number of fragments", "Anzahl Fragmente"), |s| Field::I32(&mut s.sql_max_number_fragments), Some((1.0, 100.0)), false),
    key(keys::SQL_MINIMAL_NUMBER_EQUAL_FRAGMENTS, ("Minimum number of identical frag.", "Mindestzahl gleicher Fragmente"), |s| Field::I32(&mut s.sql_minimal_number_equal_fragments), Some((0.0, 100.0)), false),
    key(keys::PEAK_DETECTOR, ("Peak detection method", "Peakerkennungsmethode"), |s| Field::Detector(&mut s.peak_detector), Some((1.0, 1000.0)), false),
];

/// Display text of a setting (English, German); panics on a key that is not in the table.
pub fn label(name: &str) -> (&'static str, &'static str) {
    KEYS.iter()
        .find(|k| k.name == name)
        .expect("a key of the table")
        .label
}

/// Bounds of a setting for an input widget; panics on a key that is not in the table.
pub fn limits(name: &str) -> (f64, f64) {
    let k = KEYS
        .iter()
        .find(|k| k.name == name)
        .expect("a key of the table");
    k.range.unwrap_or(ANY)
}

fn put<T: std::str::FromStr>(f: &mut T, v: &str) -> bool {
    v.parse().map(|x| *f = x).is_ok()
}

impl Field<'_> {
    /// Reads the app-format text.
    pub fn parse(self, v: &str) -> bool {
        match self {
            Field::F64(f) => put(f, v),
            Field::Auto(f) if v == "auto" => {
                *f = None;
                true
            }
            Field::Auto(f) => v.parse().map(|x| *f = Some(x)).is_ok(),
            Field::F32(f) => put(f, v) || variant_f32(v).map(|x| *f = x).is_some(),
            Field::I32(f) => put(f, v),
            Field::Bool(f) => put(f, v),
            Field::Ranges(f) => ranges(v).map(|r| *f = r).is_some(),
            Field::Detector(f) => detector(v).map(|d| *f = d).is_some(),
        }
    }

    /// The app-format text.
    pub fn text(self) -> String {
        match self {
            Field::F64(f) => f.to_string(),
            Field::Auto(f) => f.map_or("auto".into(), |x| x.to_string()),
            Field::F32(f) => f.to_string(),
            Field::I32(f) => f.to_string(),
            Field::Bool(f) => f.to_string(),
            Field::Ranges(f) => (f.iter().map(|(a, b)| format!("{a}:{b}")))
                .collect::<Vec<_>>()
                .join(";"),
            Field::Detector(Detector::Default) => "default".into(),
            Field::Detector(Detector::Matched { width, threshold }) => {
                format!("matched {width} {threshold}")
            }
        }
    }

    /// The old file's bytes of a `fixed` key.
    fn write(&self, out: &mut Vec<u8>) {
        match self {
            Field::F64(v) => out.extend(v.to_le_bytes()),
            // Automatic is the old default 1 there: the old version knows no other.
            Field::Auto(v) => out.extend(v.unwrap_or(1.0).to_le_bytes()),
            Field::F32(v) => out.extend(v.to_le_bytes()),
            Field::I32(v) => out.extend(v.to_le_bytes()),
            Field::Bool(v) => out.push(u8::from(**v)),
            Field::Ranges(_) | Field::Detector(_) => unreachable!("not a fixed-size key"),
        }
    }

    fn read(self, c: &mut Cursor) {
        match self {
            Field::F64(v) => *v = c.f64(),
            Field::Auto(v) => *v = Some(c.f64()),
            Field::F32(v) => *v = c.f32(),
            Field::I32(v) => *v = c.i32(),
            Field::Bool(v) => *v = c.bool(),
            Field::Ranges(_) | Field::Detector(_) => unreachable!("not a fixed-size key"),
        }
    }
}

/// `73:75;147:148`; the old app's `@Invalid()` is its empty list.
fn ranges(v: &str) -> Option<Vec<(f32, f32)>> {
    if v.is_empty() || v == "@Invalid()" {
        return Some(Vec::new());
    }
    // shortcut: a non-empty old list (a QVariantList, never seen) is not decoded and keeps the
    // default; decode it once a user's old file has one.
    (v.split(';'))
        .map(|r| {
            let (a, b) = r.split_once(':')?;
            Some((a.parse().ok()?, b.parse().ok()?))
        })
        .collect()
}

fn detector(v: &str) -> Option<Detector> {
    match v.split_whitespace().collect::<Vec<_>>()[..] {
        ["default"] => Some(Detector::Default),
        ["matched", w, t] => Some(Detector::Matched {
            width: w.parse().ok()?,
            threshold: t.parse().ok()?,
        }),
        _ => None,
    }
}

/// The old `LS_RESOLUTION=@Variant(\0\0\0\x87>\x99\x99\x9a)`: a QDataStream QVariant of type 135
/// (float), big-endian, here already unescaped to one char per byte.
fn variant_f32(v: &str) -> Option<f32> {
    let inner = v.strip_prefix("@Variant(")?.strip_suffix(')')?;
    let bytes: Vec<u8> = inner
        .chars()
        .map(|c| u8::try_from(u32::from(c)).ok())
        .collect::<Option<_>>()?;
    match bytes[..] {
        [0, 0, 0, 135, a, b, c, d] => Some(f32::from_be_bytes([a, b, c, d])),
        _ => None,
    }
}

impl GcmsSettings {
    /// Values the analysis steps can run with; the error names the first bad field.
    pub fn validate(&self) -> Result<(), String> {
        self.check().map_err(|(_, e)| e)
    }

    /// The first bad field as its table key and message. The ranges are those of [`KEYS`]; the dialog
    /// uses the same ones ([`limits`]), so it never changes a loaded value.
    fn check(&self) -> Result<(), (&'static str, String)> {
        let mut copy = self.clone();
        for k in &KEYS {
            let (lo, hi) = k.range.unwrap_or(ANY);
            let name = k.name;
            let number = |label: &str, v: f64| {
                if !v.is_finite() {
                    Err((name, format!("{label}: value out of range")))
                } else if k.range.is_some() && !(lo..=hi).contains(&v) {
                    Err((name, format!("{label}: must be between {lo} and {hi}")))
                } else if v.abs() > ANY.1 {
                    Err((name, format!("{label}: value out of range")))
                } else {
                    Ok(())
                }
            };
            match (k.field)(&mut copy) {
                Field::F64(v) => number(name, *v)?,
                Field::F32(v) => number(name, f64::from(*v))?,
                Field::I32(v) => number(name, f64::from(*v))?,
                Field::Auto(Some(v)) => number(name, *v)?,
                Field::Ranges(r) => {
                    let ok = |v: &f32| (lo..=hi).contains(&f64::from(*v));
                    if r.iter().any(|(a, b)| !ok(a) || !ok(b)) {
                        return Err((name, format!("{name}: value out of range")));
                    }
                }
                Field::Detector(Detector::Matched { width, threshold }) => {
                    number("Peak width", f64::from(*width))?;
                    // The threshold only needs to stay within the common bound.
                    if f64::from(*threshold).abs() > ANY.1 {
                        return Err((name, "Peak threshold: value out of range".into()));
                    }
                }
                Field::Auto(None) | Field::Bool(_) | Field::Detector(Detector::Default) => {}
            }
        }
        Ok(())
    }

    /// Sets one setting from its app-format text. `None`: no such key; `Some(false)`: unreadable text.
    pub fn set(&mut self, key: &str, text: &str) -> Option<bool> {
        let k = KEYS.iter().find(|k| k.name == key)?;
        Some((k.field)(self).parse(text))
    }

    /// Every setting as `key`, app-format text, in table order.
    pub fn lines(&self) -> Vec<(&'static str, String)> {
        let mut copy = self.clone();
        (KEYS.iter())
            .map(|k| (k.name, (k.field)(&mut copy).text()))
            .collect()
    }

    /// Loaded values that fail [`Self::validate`] go back to their default, one field at a time
    /// (each message is pushed to `bad`); all of them if a message names no field.
    pub fn sanitized(mut self, bad: &mut Vec<String>) -> Self {
        // Each round resets one field, so this ends; the bound only guards against a field whose
        // default does not pass.
        for _ in 0..=KEYS.len() {
            let Err((name, e)) = self.check() else {
                return self;
            };
            bad.push(e);
            let Some(k) = KEYS.iter().find(|k| k.name == name) else {
                break;
            };
            let default = (k.field)(&mut Self::default()).text();
            (k.field)(&mut self).parse(&default);
        }
        Self::default()
    }

    /// Old `toStream` (§5.2): 191 + 8n bytes, no header.
    pub fn to_bytes(&self) -> Vec<u8> {
        let mut out = Vec::with_capacity(FIXED_LEN + ION_LEN * self.quant_excl_ions.len());
        let mut copy = self.clone();
        for k in KEYS.iter().filter(|k| k.fixed) {
            (k.field)(&mut copy).write(&mut out);
        }
        out.extend((self.quant_excl_ions.len() as u32).to_le_bytes());
        for (a, b) in &self.quant_excl_ions {
            out.extend(a.to_le_bytes());
            out.extend(b.to_le_bytes());
        }
        out.extend(self.checksum().to_le_bytes());
        out
    }

    fn checksum(&self) -> f64 {
        self.quant_excl_ions.len() as f64
            + self.sf_peak_threshold_begin
            + self.ls_pure_factor
            + self.an_peak_threshold_begin
    }

    /// Reads an old settings file (raw, or gzip-wrapped as algorithms.md §1 suggests) and returns
    /// `self` with the stored values replaced; the values the file does not hold stay as they are.
    /// Rejects wrong sizes, a wrong checksum and values the analysis cannot use.
    pub fn from_bytes(&self, bytes: &[u8]) -> Result<Self, String> {
        let mut unpacked = Vec::new();
        // A raw file can start with these bytes by chance (they are the low bytes of the first
        // value), so a failed gzip read of a file of raw size falls back to the raw layout.
        let gz = bytes.starts_with(&[0x1f, 0x8b, 0x08]);
        let bytes = match gz.then(|| {
            GzDecoder::new(bytes)
                .take(MAX_LEN + 1)
                .read_to_end(&mut unpacked)
        }) {
            Some(Ok(_)) if unpacked.len() as u64 > MAX_LEN => {
                return Err("settings file too large".into());
            }
            Some(Ok(_)) => &unpacked[..],
            Some(Err(e)) if bytes.len() < FIXED_LEN => {
                return Err(format!("not a readable gzip file: {e}"));
            }
            _ => bytes,
        };
        if bytes.len() < FIXED_LEN {
            return Err(format!("settings file too short ({} bytes)", bytes.len()));
        }
        let n = u32::from_le_bytes(bytes[COUNT_AT..COUNT_AT + 4].try_into().expect("4 bytes"));
        let expect = (n as usize)
            .checked_mul(ION_LEN)
            .and_then(|b| b.checked_add(FIXED_LEN));
        if expect != Some(bytes.len()) {
            return Err(format!(
                "settings file has {} bytes, not 191 + 8 x {n}",
                bytes.len()
            ));
        }
        let mut c = Cursor(bytes);
        let mut s = self.clone();
        for k in KEYS.iter().filter(|k| k.fixed) {
            (k.field)(&mut s).read(&mut c);
        }
        c.take::<4>();
        s.quant_excl_ions = (0..n).map(|_| (c.f32(), c.f32())).collect();
        let check = c.f64();
        if (check - s.checksum())
            .abs()
            .partial_cmp(&(1e-9 * (1.0 + check.abs())))
            != Some(std::cmp::Ordering::Less)
        {
            return Err("settings file checksum does not match".into());
        }
        s.validate()?;
        Ok(s)
    }

    /// The settings of every analysis step, the one derivation from these settings. `AN_NUMBER_BINS`
    /// is read once here, so the peak and the deconvolution step cannot differ.
    pub fn analysis(&self) -> AnalysisSettings {
        let number_bins = self.an_number_bins.max(1) as u32;
        AnalysisSettings {
            peaks: PeakSettings {
                threshold_begin: self.an_peak_threshold_begin as f32,
                threshold_end: self.an_peak_threshold_end as f32,
                min_peak_height: self.an_min_peak_height as f32,
                number_bins,
                detector: self.peak_detector,
            },
            deconv: DeconvSettings {
                width: self.an_deconvolution_width,
                number_bins,
                min_peak_number: self.an_min_peak_number.max(1) as usize,
                min_height_bp: self.an_min_height_bp,
                modeled_integration: self.an_modeled_values_for_integration,
                excluded_model_ions: self.quant_excl_ions.clone(),
            },
            quant: QuantSettings {
                ion_number: self.quant_ion_number.max(1) as usize,
                min_dist: self.quant_min_dist_ions as f32,
                quality_cutoff: self.quant_quality_cutoff as f32,
                excluded: self.quant_excl_ions.clone(),
            },
            // The search cuts off at `LS_CUTOFF` everywhere, a deliberate deviation from 3.3 (algorithms.md §5.8).
            search: SearchSettings {
                ri_diff: self.ls_ri_diff,
                pure_factor: self.ls_pure_factor,
                impure_factor: self.ls_impure_factor,
                scale_ion: self.ls_scale_factor_ion,
                scale_value: self.ls_scale_factor_value,
                no_match: self.ls_scale_factor_no_match,
                use_ri: self.ls_use_ri,
                use_spec: self.ls_use_spec,
                mass_filter: self.ls_mass_filter.clone(),
                cutoff: self.ls_cutoff,
                max_fragments: self.sql_max_number_fragments.max(1) as usize,
                min_equal_fragments: self.sql_minimal_number_equal_fragments.max(0) as usize,
            },
            fwhm: self.sf_fwhm,
        }
    }

    /// Settings of a batch scan: the library is prepared with the search settings, which cut hits off at
    /// `LS_CUTOFF` like the old `BatchQuantificator`.
    pub fn scan(&self, calibration: Option<RiCalibration>) -> ScanSettings {
        ScanSettings {
            analysis: self.analysis(),
            calibration,
            shift: None,
            untargeted: false,
            save_cmp: None,
            baseline_correction: false,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The defaults of file-formats.md §5.1 and the old defaults of each step.
    #[test]
    fn defaults_match_old_static_values() {
        let s = GcmsSettings::default();
        assert_eq!(s.analysis().peaks, PeakSettings::default());
        // SF_FWHM reaches every centroiding (#293).
        let f = GcmsSettings {
            sf_fwhm: 0.3,
            ..s.clone()
        };
        let centroid = crate::centroid::CentroidSettings::from(&f);
        assert_eq!((f.analysis().fwhm, centroid.fwhm), (0.3, 0.3));
        // Except the deconvolution width: automatic (#279), written to the old file as its 1.
        assert_eq!(
            s.analysis().deconv,
            DeconvSettings {
                width: None,
                ..Default::default()
            }
        );
        let back = GcmsSettings::default().from_bytes(&s.to_bytes()).unwrap();
        assert_eq!(back.an_deconvolution_width, Some(1.0));
        assert_eq!(s.analysis().search, SearchSettings::default());
        assert_eq!(
            (
                s.sf_peak_threshold_end,
                s.sf_fwhm,
                s.sf_centroid_max_base_dev
            ),
            (-10.0, 0.5, 30.0)
        );
        assert_eq!((s.ls_resolution, s.ls_use_scaled_lib), (0.3, true));
        assert_eq!(s.ls_cutoff, SearchSettings::default().cutoff);
        assert_eq!(
            (
                s.quant_ion_number,
                s.quant_min_dist_ions,
                s.quant_quality_cutoff
            ),
            (3, 0.5, 1.0)
        );
        s.validate().unwrap();
        assert_eq!(
            (
                s.analysis().search.max_fragments,
                s.analysis().search.min_equal_fragments
            ),
            (5, 3)
        );
        // Not in the old file, so only the dialog can set them; it keeps them in range.
        for (max, min) in [(0, 3), (101, 3), (5, -1), (5, 101)] {
            let bad = GcmsSettings {
                sql_max_number_fragments: max,
                sql_minimal_number_equal_fragments: min,
                ..Default::default()
            };
            assert!(
                bad.validate().unwrap_err().starts_with("SQL_"),
                "{max} {min}"
            );
        }
    }

    /// A value that differs from the default and passes `validate`.
    fn change(k: &Key, s: &mut GcmsSettings) {
        match (k.field)(s) {
            Field::F64(f) => *f = k.range.map_or(3.5, |(lo, hi)| lo + (hi - lo) / 3.0),
            Field::Auto(f) => *f = Some(2.5),
            Field::F32(f) => *f = 0.25,
            Field::I32(f) => *f = k.range.map_or(7, |(lo, _)| lo as i32 + 1),
            Field::Bool(f) => *f = !*f,
            Field::Ranges(f) => *f = vec![(73.0, 75.0), (147.0, 148.5)],
            Field::Detector(f) => {
                *f = Detector::Matched {
                    width: 12,
                    threshold: 4,
                }
            }
        }
    }

    /// Puts `v` where the key's range applies: the number, the matched width, or both ends of a range.
    fn put_bound(k: &Key, s: &mut GcmsSettings, v: f64) {
        match (k.field)(s) {
            Field::F64(f) => *f = v,
            Field::F32(f) => *f = v as f32,
            Field::I32(f) => *f = v as i32,
            Field::Auto(f) => *f = Some(v),
            Field::Detector(d) => {
                *d = Detector::Matched {
                    width: v as u32,
                    threshold: 0,
                }
            }
            Field::Ranges(r) => *r = vec![(v as f32, v as f32)],
            Field::Bool(_) => unreachable!("a flag has no range"),
        }
    }

    /// Each key round-trips through the app format and, if it has a slot there, the old file.
    #[test]
    fn every_key_round_trips() {
        // The old file holds the width, so a loaded one is set.
        let base = GcmsSettings {
            an_deconvolution_width: Some(1.0),
            ..Default::default()
        };
        let names: std::collections::HashSet<_> = KEYS.iter().map(|k| k.name).collect();
        assert_eq!(names.len(), KEYS.len(), "a key twice");
        assert!(
            KEYS.iter()
                .all(|k| !k.label.0.is_empty() && !k.label.1.is_empty())
        );
        for k in &KEYS {
            // Default: the text of every key reads back to the same value.
            let mut back = GcmsSettings::default();
            for (name, text) in GcmsSettings::default().lines() {
                assert_eq!(back.set(name, &text), Some(true), "{name}");
            }
            assert_eq!(back, GcmsSettings::default(), "{}", k.name);

            let mut s = base.clone();
            change(k, &mut s);
            assert_ne!(s, base, "{}", k.name);
            s.validate().unwrap_or_else(|e| panic!("{}: {e}", k.name));
            let mut app = base.clone();
            for (name, text) in s.lines() {
                assert_eq!(app.set(name, &text), Some(true), "{name}");
            }
            assert_eq!(app, s, "app format {}", k.name);
            let old = base.from_bytes(&s.to_bytes()).unwrap();
            if k.fixed {
                assert_eq!(old, s, "old file {}", k.name);
            } else if k.name != "QUANT_EXCL_IONS" {
                assert_eq!(old, base, "{} is not in the old file", k.name);
            }
        }
    }

    #[test]
    fn set_names_unknown_keys_and_bad_text() {
        assert_eq!(GcmsSettings::default().set("NO_SUCH_KEY", "1"), None);
        assert_eq!(
            GcmsSettings::default().set("LS_RI_DIFF", "abc"),
            Some(false)
        );
    }

    /// The bounds of the table are the bounds of `validate`, at the edge and one step beyond.
    #[test]
    fn every_range_is_enforced() {
        for k in KEYS.iter().filter(|k| k.range.is_some()) {
            let (lo, hi) = k.range.unwrap();
            for v in [lo, hi] {
                let mut s = GcmsSettings::default();
                put_bound(k, &mut s, v);
                s.validate()
                    .unwrap_or_else(|e| panic!("{} {v}: {e}", k.name));
            }
            for v in [lo - 1.0, hi + 1.0] {
                let mut s = GcmsSettings::default();
                put_bound(k, &mut s, v);
                assert!(
                    s.check().is_err_and(|(name, _)| name == k.name),
                    "{} {v}",
                    k.name
                );
            }
        }
        // A bad value is reset alone, with its message.
        let mut bad = Vec::new();
        let s = GcmsSettings {
            an_number_bins: 0,
            ls_pure_factor: 0.7,
            peak_detector: Detector::Matched {
                width: 0,
                threshold: 5,
            },
            ..Default::default()
        };
        let s = s.sanitized(&mut bad);
        assert_eq!(
            bad,
            [
                "AN_NUMBER_BINS: must be between 1 and 100",
                "Peak width: must be between 1 and 1000"
            ]
        );
        assert_eq!(
            s,
            GcmsSettings {
                ls_pure_factor: 0.7,
                ..Default::default()
            }
        );
    }

    /// `peak` and `deconv` read the apex steps from one place, and the excluded ions reach the model peak
    /// and the quantification.
    #[test]
    fn shared_values_reach_every_step() {
        let s = GcmsSettings {
            an_number_bins: 7,
            quant_excl_ions: vec![(73.0, 73.0)],
            ..Default::default()
        };
        assert_eq!(
            (
                s.analysis().peaks.number_bins,
                s.analysis().deconv.number_bins
            ),
            (7, 7)
        );
        assert_eq!(
            s.analysis().deconv.excluded_model_ions,
            s.analysis().quant.excluded
        );
    }

    #[test]
    fn scan_and_quant_settings_follow_the_fields() {
        let s = GcmsSettings {
            ls_cutoff: 0.6,
            quant_ion_number: 0,
            quant_min_dist_ions: 1.5,
            quant_quality_cutoff: 2.5,
            quant_excl_ions: vec![(73.0, 73.0)],
            ..Default::default()
        };
        let q = s.analysis().quant;
        assert_eq!((q.ion_number, q.min_dist, q.quality_cutoff), (1, 1.5, 2.5));
        assert_eq!(q.excluded, [(73.0, 73.0)]);
        let cal = RiCalibration::new([(0.0, 0.0), (60.0, 100.0)]).unwrap();
        let scan = s.scan(Some(cal.clone()));
        assert_eq!(
            (scan.analysis.peaks, scan.calibration),
            (s.analysis().peaks, Some(cal))
        );
        // One cutoff for the library search and the scan (#399).
        assert_eq!(scan.analysis.search, s.analysis().search);
        assert_eq!(scan.analysis.search.cutoff, s.ls_cutoff);
    }

    /// §5.2: 191 bytes for the defaults, checksum 20.5.
    #[test]
    fn default_file_is_191_bytes_with_checksum_20_5() {
        let b = GcmsSettings::default().to_bytes();
        assert_eq!(b.len(), 191);
        assert_eq!(f64::from_le_bytes(b[183..].try_into().unwrap()), 20.5);
        assert_eq!(f64::from_le_bytes(b[0..8].try_into().unwrap()), 10.0);
        assert_eq!(f32::from_le_bytes(b[152..156].try_into().unwrap()), 0.3);
    }

    #[test]
    fn round_trips_with_excluded_ions_and_gzip() {
        let s = GcmsSettings {
            an_peak_threshold_begin: 12.5,
            an_number_bins: 7,
            an_deconvolution_width: Some(2.5),
            ls_use_ri: false,
            quant_excl_ions: vec![(73.0, 73.0), (147.0, 149.5)],
            ..Default::default()
        };
        let b = s.to_bytes();
        assert_eq!(b.len(), 191 + 16);
        assert_eq!(GcmsSettings::default().from_bytes(&b).unwrap(), s);

        use std::io::Write;
        let mut z = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::default());
        z.write_all(&b).unwrap();
        assert_eq!(
            GcmsSettings::default()
                .from_bytes(&z.finish().unwrap())
                .unwrap(),
            s
        );
    }

    #[test]
    fn keeps_values_the_file_does_not_hold() {
        let base = GcmsSettings {
            an_min_height_bp: 0.25,
            ls_mass_filter: vec![(1.0, 2.0)],
            ls_cutoff: 0.85,
            sql_max_number_fragments: 7,
            sql_minimal_number_equal_fragments: 2,
            // The file holds the width, so a loaded one is set.
            an_deconvolution_width: Some(1.0),
            ..Default::default()
        };
        let got = base
            .from_bytes(&GcmsSettings::default().to_bytes())
            .unwrap();
        assert_eq!(got, base);
    }

    #[test]
    fn rejects_bad_files() {
        let good = GcmsSettings::default().to_bytes();
        let d = GcmsSettings::default();
        assert!(d.from_bytes(&[]).is_err());
        assert!(d.from_bytes(&good[..190]).is_err());
        let mut long = good.clone();
        long.push(0);
        assert!(d.from_bytes(&long).is_err());
        // Huge excluded-ion count with a short file: no allocation, no panic.
        let mut huge = good.clone();
        huge[179..183].copy_from_slice(&u32::MAX.to_le_bytes());
        assert!(d.from_bytes(&huge).is_err());
        let mut sum = good.clone();
        sum[190] ^= 0x40;
        assert!(d.from_bytes(&sum).unwrap_err().contains("checksum"));
        // Valid checksum but NaN and out-of-range values.
        for edit in [
            |s: &mut GcmsSettings| s.an_peak_threshold_end = f64::NAN,
            |s: &mut GcmsSettings| s.an_number_bins = 0,
            |s: &mut GcmsSettings| s.an_number_bins = i32::MAX,
            |s: &mut GcmsSettings| s.an_deconvolution_width = Some(0.0),
            |s: &mut GcmsSettings| s.an_min_peak_number = -1,
            |s: &mut GcmsSettings| s.ls_ri_diff = 0.0,
            |s: &mut GcmsSettings| s.quant_excl_ions = vec![(f32::NAN, 1.0)],
            |s: &mut GcmsSettings| s.quant_excl_ions = vec![(-1.0, 73.0)],
        ] {
            let mut s = GcmsSettings::default();
            edit(&mut s);
            assert!(d.from_bytes(&s.to_bytes()).is_err(), "{s:?}");
        }
        // Truncated and oversized gzip.
        assert!(d.from_bytes(&[0x1f, 0x8b, 8, 0]).is_err());
        let mut z = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::default());
        std::io::Write::write_all(&mut z, &vec![0u8; MAX_LEN as usize + 1]).unwrap();
        assert!(
            d.from_bytes(&z.finish().unwrap())
                .unwrap_err()
                .contains("too large")
        );
    }

    /// A raw file whose first value starts with the gzip magic bytes still loads.
    #[test]
    fn raw_file_starting_with_gzip_magic_loads() {
        let s = GcmsSettings {
            sf_peak_threshold_begin: f64::from_le_bytes([0x1f, 0x8b, 0x08, 0, 0, 0, 0x24, 0x40]),
            an_deconvolution_width: Some(1.0),
            ..Default::default()
        };
        let b = s.to_bytes();
        assert_eq!(b[..3], [0x1f, 0x8b, 0x08]);
        assert_eq!(GcmsSettings::default().from_bytes(&b).unwrap(), s);
    }
}
