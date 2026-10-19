//! The old `GCMSSettings` (file-formats.md §5): all values in one struct, the old file format, and the
//! conversion to the settings of each analysis step. Files are untrusted: `from_bytes` checks size,
//! checksum and every value before anything is used.

use std::io::Read;

use flate2::read::GzDecoder;

use crate::deconv::DeconvSettings;
use crate::peaks::PeakSettings;
use crate::search::SearchSettings;

/// Old `GCMSSettings`. Defaults are the old static values (file-formats.md §5.1).
/// Fields marked "stored only" are kept for the file but no step of the rewrite reads them yet.
#[derive(Debug, Clone, PartialEq)]
pub struct GcmsSettings {
    // Spectrum filter (stored only).
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
    pub an_deconvolution_width: f64,
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
            an_deconvolution_width: 1.0,
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

impl GcmsSettings {
    /// Values the analysis steps can run with; the error names the first bad field.
    pub fn validate(&self) -> Result<(), String> {
        let finite = |name: &str, v: f64| {
            if v.is_finite() && v.abs() <= 1e9 {
                Ok(())
            } else {
                Err(format!("{name}: value out of range"))
            }
        };
        for (name, v) in [
            ("SF_PEAK_THRESHOLD_BEGIN", self.sf_peak_threshold_begin),
            ("SF_PEAK_THRESHOLD_END", self.sf_peak_threshold_end),
            ("SF_MAX_BASE_DEV", self.sf_max_base_dev),
            ("SF_FWHM", self.sf_fwhm),
            (
                "SF_CENTROID_THRESHOLD_BEGIN",
                self.sf_centroid_threshold_begin,
            ),
            ("SF_CENTROID_THRESHOLD_END", self.sf_centroid_threshold_end),
            ("SF_CENTROID_MAX_BASE_DEV", self.sf_centroid_max_base_dev),
            ("AN_PEAK_THRESHOLD_BEGIN", self.an_peak_threshold_begin),
            ("AN_PEAK_THRESHOLD_END", self.an_peak_threshold_end),
            ("AN_MAX_BASE_DEV", self.an_max_base_dev),
            ("AN_MIN_PEAK_HEIGHT", self.an_min_peak_height),
            ("LS_PURE_FACTOR", self.ls_pure_factor),
            ("LS_IMPURE_FACTOR", self.ls_impure_factor),
            ("LS_SCALE_FACTOR_ION", self.ls_scale_factor_ion),
            ("LS_SCALE_FACTOR_VALUE", self.ls_scale_factor_value),
            ("LS_SCALE_FACTOR_NO_MATCH", self.ls_scale_factor_no_match),
            ("QUANT_MIN_DIST_IONS", self.quant_min_dist_ions),
            ("QUANT_QUALITY_CUTOFF", self.quant_quality_cutoff),
            ("LS_RESOLUTION", self.ls_resolution as f64),
        ] {
            finite(name, v)?;
        }
        let ranged = |name: &str, v: f64, lo: f64, hi: f64| {
            if (lo..=hi).contains(&v) {
                Ok(())
            } else {
                Err(format!("{name}: must be between {lo} and {hi}"))
            }
        };
        ranged("AN_MIN_PEAK_HEIGHT", self.an_min_peak_height, 0.0, 1e9)?;
        ranged(
            "AN_DECONVOLUTION_WIDTH",
            self.an_deconvolution_width,
            0.1,
            100.0,
        )?;
        ranged("AN_NUMBER_BINS", self.an_number_bins as f64, 1.0, 100.0)?;
        ranged(
            "AN_MIN_PEAK_NUMBER",
            self.an_min_peak_number as f64,
            1.0,
            10_000.0,
        )?;
        ranged("AN_MIN_HEIGHT_BP", self.an_min_height_bp, 0.0, 1.0)?;
        ranged("LS_RI_DIFF", self.ls_ri_diff, 1e-6, 1e9)?;
        ranged(
            "QUANT_ION_NUMBER",
            self.quant_ion_number as f64,
            0.0,
            10_000.0,
        )?;
        for (name, v) in [
            ("LS_PURE_FACTOR", self.ls_pure_factor),
            ("LS_IMPURE_FACTOR", self.ls_impure_factor),
        ] {
            ranged(name, v, 0.0, 1e9)?;
        }
        for (name, ranges) in [
            ("QUANT_EXCL_IONS", &self.quant_excl_ions),
            ("LS_MASS_FILTER", &self.ls_mass_filter),
        ] {
            if ranges.iter().any(|(a, b)| !a.is_finite() || !b.is_finite()) {
                return Err(format!("{name}: value out of range"));
            }
        }
        Ok(())
    }

    /// Old `toStream` (§5.2): 191 + 8n bytes, no header.
    pub fn to_bytes(&self) -> Vec<u8> {
        let mut out = Vec::with_capacity(FIXED_LEN + ION_LEN * self.quant_excl_ions.len());
        for v in [
            self.sf_peak_threshold_begin,
            self.sf_peak_threshold_end,
            self.sf_max_base_dev,
            self.sf_fwhm,
            self.sf_centroid_threshold_begin,
            self.sf_centroid_threshold_end,
            self.sf_centroid_max_base_dev,
            self.an_peak_threshold_begin,
            self.an_peak_threshold_end,
            self.an_max_base_dev,
            self.an_min_peak_height,
            self.an_deconvolution_width,
        ] {
            out.extend(v.to_le_bytes());
        }
        out.extend(self.an_number_bins.to_le_bytes());
        out.extend(self.an_min_peak_number.to_le_bytes());
        for v in [
            self.ls_ri_diff,
            self.ls_pure_factor,
            self.ls_impure_factor,
            self.ls_scale_factor_ion,
            self.ls_scale_factor_value,
            self.ls_scale_factor_no_match,
        ] {
            out.extend(v.to_le_bytes());
        }
        out.extend(self.ls_resolution.to_le_bytes());
        out.extend([self.ls_use_ri, self.ls_use_scaled_lib, self.ls_use_spec].map(u8::from));
        out.extend(self.quant_ion_number.to_le_bytes());
        out.extend(self.quant_min_dist_ions.to_le_bytes());
        out.extend(self.quant_quality_cutoff.to_le_bytes());
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
        let bytes = if bytes.starts_with(&[0x1f, 0x8b]) {
            GzDecoder::new(bytes)
                .take(MAX_LEN + 1)
                .read_to_end(&mut unpacked)
                .map_err(|e| format!("not a readable gzip file: {e}"))?;
            if unpacked.len() as u64 > MAX_LEN {
                return Err("settings file too large".into());
            }
            &unpacked[..]
        } else {
            bytes
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
        let mut s = Self {
            sf_peak_threshold_begin: c.f64(),
            sf_peak_threshold_end: c.f64(),
            sf_max_base_dev: c.f64(),
            sf_fwhm: c.f64(),
            sf_centroid_threshold_begin: c.f64(),
            sf_centroid_threshold_end: c.f64(),
            sf_centroid_max_base_dev: c.f64(),
            an_peak_threshold_begin: c.f64(),
            an_peak_threshold_end: c.f64(),
            an_max_base_dev: c.f64(),
            an_min_peak_height: c.f64(),
            an_deconvolution_width: c.f64(),
            an_number_bins: c.i32(),
            an_min_peak_number: c.i32(),
            ls_ri_diff: c.f64(),
            ls_pure_factor: c.f64(),
            ls_impure_factor: c.f64(),
            ls_scale_factor_ion: c.f64(),
            ls_scale_factor_value: c.f64(),
            ls_scale_factor_no_match: c.f64(),
            ls_resolution: c.f32(),
            ls_use_ri: c.bool(),
            ls_use_scaled_lib: c.bool(),
            ls_use_spec: c.bool(),
            quant_ion_number: c.i32(),
            quant_min_dist_ions: c.f64(),
            quant_quality_cutoff: c.f64(),
            ..self.clone()
        };
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

    pub fn peak(&self) -> PeakSettings {
        PeakSettings {
            threshold_begin: self.an_peak_threshold_begin as f32,
            threshold_end: self.an_peak_threshold_end as f32,
            min_peak_height: self.an_min_peak_height as f32,
            number_bins: self.an_number_bins.max(1) as u32,
        }
    }

    pub fn deconv(&self) -> DeconvSettings {
        DeconvSettings {
            width: self.an_deconvolution_width,
            number_bins: self.an_number_bins.max(1) as u32,
            min_peak_number: self.an_min_peak_number.max(1) as usize,
            min_height_bp: self.an_min_height_bp,
            modeled_integration: self.an_modeled_values_for_integration,
            excluded_model_ions: self.quant_excl_ions.clone(),
        }
    }

    /// The GUI searches with cutoff 0 and colours by score (algorithms.md §5.8), so `LS_CUTOFF` is not used.
    pub fn search(&self) -> SearchSettings {
        SearchSettings {
            ri_diff: self.ls_ri_diff,
            pure_factor: self.ls_pure_factor,
            impure_factor: self.ls_impure_factor,
            scale_ion: self.ls_scale_factor_ion,
            scale_value: self.ls_scale_factor_value,
            no_match: self.ls_scale_factor_no_match,
            use_ri: self.ls_use_ri,
            use_spec: self.ls_use_spec,
            mass_filter: self.ls_mass_filter.clone(),
            cutoff: 0.0,
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
        assert_eq!(s.peak(), PeakSettings::default());
        assert_eq!(s.deconv(), DeconvSettings::default());
        assert_eq!(
            s.search(),
            SearchSettings {
                cutoff: 0.0,
                ..Default::default()
            }
        );
        assert_eq!(
            (
                s.sf_peak_threshold_end,
                s.sf_fwhm,
                s.sf_centroid_max_base_dev
            ),
            (-10.0, 0.5, 30.0)
        );
        assert_eq!((s.ls_resolution, s.ls_use_scaled_lib), (0.3, true));
        assert_eq!(
            (
                s.quant_ion_number,
                s.quant_min_dist_ions,
                s.quant_quality_cutoff
            ),
            (3, 0.5, 1.0)
        );
        s.validate().unwrap();
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
            |s: &mut GcmsSettings| s.an_deconvolution_width = 0.0,
            |s: &mut GcmsSettings| s.an_min_peak_number = -1,
            |s: &mut GcmsSettings| s.ls_ri_diff = 0.0,
            |s: &mut GcmsSettings| s.quant_excl_ions = vec![(f32::NAN, 1.0)],
        ] {
            let mut s = GcmsSettings::default();
            edit(&mut s);
            assert!(d.from_bytes(&s.to_bytes()).is_err(), "{s:?}");
        }
        // Truncated or oversized gzip.
        assert!(d.from_bytes(&[0x1f, 0x8b, 8, 0]).is_err());
    }
}
