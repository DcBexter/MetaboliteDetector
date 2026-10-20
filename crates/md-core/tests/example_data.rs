//! Reads every real run in `example-data/` (fetched by `example-data/fetch.sh` or `fetch.ps1`).
//! Passes without checking anything when the data has not been downloaded.

use std::path::{Path, PathBuf};
use std::time::Instant;

use md_core::deconv::{DeconvSettings, deconvolute};
use md_core::peaks::{IonTraces, PeakSettings};

fn runs() -> Vec<PathBuf> {
    let dir = Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../../example-data/TargetSearchData/inst/gc-ms-data");
    let Ok(entries) = std::fs::read_dir(&dir) else {
        eprintln!(
            "skipped: {} not found, run example-data/fetch.*",
            dir.display()
        );
        return Vec::new();
    };
    let runs: Vec<_> = entries
        .map(|e| e.unwrap().path())
        .filter(|p| p.extension().is_some_and(|e| e.eq_ignore_ascii_case("cdf")))
        .collect();
    assert!(!runs.is_empty(), "no .cdf files in {}", dir.display());
    runs
}

fn name(path: &Path) -> std::borrow::Cow<'_, str> {
    path.file_name().unwrap().to_string_lossy()
}

#[test]
fn reads_all_example_runs() {
    for path in runs() {
        let t = Instant::now();
        let c = md_core::andi::read(&path).unwrap_or_else(|e| panic!("{}: {e}", path.display()));
        let tic = c.tic();
        assert!(c.scan_count() > 100, "{}", path.display());
        assert!(tic.iter().any(|&v| v > 0.0), "{}", path.display());
        let mz = |s| c.spectrum(s).0.to_vec();
        assert!((0..c.scan_count()).all(|s| mz(s).is_sorted()));
        eprintln!(
            "{}: {} scans, {} points, {:.1}-{:.1} min, read+TIC {:?}",
            name(&path),
            c.scan_count(),
            c.mz.len(),
            c.scan_times_s[0] / 60.0,
            c.scan_times_s[c.scan_count() - 1] / 60.0,
            t.elapsed()
        );
    }
}

/// The old app's `.bin`/`.idx` cache of each run (made by
/// docs/reverse-engineering/tools/diskscan/make-cache.sh into `$MD_OLD_CACHE`) holds the same spectra.
#[test]
fn old_cache_matches_cdf() {
    let Some(dir) = std::env::var_os("MD_OLD_CACHE").map(PathBuf::from) else {
        eprintln!("skipped: MD_OLD_CACHE not set");
        return;
    };
    for path in runs() {
        let bin = dir.join(path.file_name().unwrap()).with_extension("bin");
        let old =
            md_core::diskscan::read(&bin).unwrap_or_else(|e| panic!("{}: {e}", bin.display()));
        let cdf = md_core::andi::read(&path).unwrap();
        let n = name(&path);
        assert_eq!(
            (&old.offsets, &old.mz, &old.intensity),
            (&cdf.offsets, &cdf.mz, &cdf.intensity),
            "{n}"
        );
        // The old importer kept times as f32 seconds before scaling them to ms.
        let mut times = old.scan_times_s.iter().zip(&cdf.scan_times_s);
        assert!(times.all(|(a, b)| (a - b).abs() <= 1e-7 * b.abs()), "{n}");
        // Export to NetCDF and read back.
        let out = std::env::temp_dir().join(format!("md-core-export-{}-{n}", std::process::id()));
        md_core::andi::write(&out, &old).unwrap();
        let back = md_core::andi::read(&out);
        std::fs::remove_file(&out).ok();
        assert_eq!(back.unwrap(), old, "{n}");
    }
}

#[test]
fn detects_peaks_in_example_runs() {
    for path in runs() {
        let c = md_core::andi::read(&path).unwrap();
        let t = Instant::now();
        let mut traces = IonTraces::new(&c).unwrap();
        let index = t.elapsed();
        let peaks = traces.detect_peaks(&PeakSettings::default());
        let total = t.elapsed();
        let per_ion = peaks.len() as f64 / traces.len() as f64;
        eprintln!(
            "{}: {} traces, {} peaks ({per_ion:.1} per trace), index {index:?}, index+peaks {total:?}",
            name(&path),
            traces.len(),
            peaks.len(),
        );
        // About 2 per trace on these runs (all clear peaks found, see docs/images); far outside
        // this range means a broken detector, not different data.
        assert!((1.0..200.0).contains(&per_ion), "{}", path.display());
        for p in &peaks {
            assert!(
                p.first_scan < p.apex_scan && p.apex_scan < p.last_scan,
                "{p:?}"
            );
            assert!(p.signal_to_noise >= 20.0, "{p:?}");
        }
    }
}

#[test]
fn deconvolutes_example_runs() {
    for path in runs() {
        let c = md_core::andi::read(&path).unwrap();
        let t = Instant::now();
        let mut traces = IonTraces::new(&c).unwrap();
        let peaks = traces.detect_peaks(&PeakSettings::default());
        let s = DeconvSettings::default();
        let comps = deconvolute(&traces, &peaks, &c.scan_times_s, &s);
        let total = t.elapsed();
        let ions: usize = comps.iter().map(|k| k.spectrum.len()).sum();
        eprintln!(
            "{}: {} peaks -> {} components ({:.1} ions each), index+peaks+deconvolution {total:?}",
            name(&path),
            peaks.len(),
            comps.len(),
            ions as f64 / comps.len().max(1) as f64,
        );
        // A derivatised plant extract has dozens of substances; none or thousands means a broken step.
        assert!((10..500).contains(&comps.len()), "{}", path.display());
        for k in &comps {
            // Like the old code, duplicated peaks count towards the minimum, distinct masses may be fewer.
            assert!(
                k.peaks.len() >= s.min_peak_number && !k.spectrum.is_empty(),
                "{k:?}"
            );
            assert!(
                k.first_scan <= k.apex_scan && k.apex_scan <= k.last_scan,
                "{k:?}"
            );
        }
    }
}

/// Tab-separated rows of a TargetSearchData text file, header skipped.
const DENSITY: usize = 40;
const THRESHOLD: i32 = 100;

fn table(path: &Path) -> Vec<Vec<String>> {
    std::fs::read_to_string(path)
        .unwrap()
        .lines()
        .skip(1)
        .map(|l| l.split('\t').map(str::to_owned).collect())
        .collect()
}

#[test]
fn ri_of_markers_matches_target_search() {
    let runs = runs();
    let Some(dir) = runs.first().and_then(|p| p.parent()) else {
        return;
    };
    // FAME markers: time window (s) and RI; TargetSearch finds them on m/z 87. The RI files list
    // TargetSearch's RI for each scan, so the comparison covers interpolation and extrapolation.
    let windows: Vec<(f64, f64, f64)> = table(&dir.join("rimLimits.txt"))
        .iter()
        .map(|r| {
            (
                r[0].parse().unwrap(),
                r[1].parse().unwrap(),
                r[2].parse().unwrap(),
            )
        })
        .collect();
    for path in &runs {
        let c = md_core::andi::read(path).unwrap();
        let peaks = md_core::ri::marker_peaks(&c, 87.0, DENSITY, THRESHOLD);
        let points = md_core::ri::window_points(&peaks, &windows);
        assert_eq!(points.len(), windows.len(), "{}: {peaks:?}", name(path));
        let cal = md_core::ri::RiCalibration::new(points.iter().copied()).unwrap();
        let reference = table(&dir.join(format!(
            "RI_{}.txt",
            path.file_stem().unwrap().to_string_lossy()
        )));
        let worst = reference
            .iter()
            .map(|r| {
                let (rt, ri): (f64, f64) = (r[0].parse().unwrap(), r[2].parse().unwrap());
                (cal.ri(rt) - ri).abs()
            })
            .fold(0.0, f64::max);
        eprintln!(
            "{}: markers {points:?}, worst RI difference {worst:.1}",
            name(path)
        );
        // TargetSearch put the middle marker of 7235eg06 at 321.1 s, where m/z 87 reaches 36 512;
        // the marker peak at 311.1 s (as in every other run) reaches 221 839.
        if name(path) == "7235eg06.cdf" {
            assert!((points[1].0 - 311.109).abs() < 1e-6, "{points:?}");
        } else {
            assert!(worst < 0.01, "{}: RI off by {worst}", name(path));
        }
    }
}

/// RI search with the sample library `lib_aa.lbr` (alkane RI; `MD_LBR_SAMPLE` as in `lbr_sample.rs`).
/// The example runs have FAME markers but no alkanes, so their alkane calibration is bridged over two
/// amino acids: the FAME RI in TargetSearch's `library.txt` gives the time in the run, `lib_aa` the RI.
/// Amino acids between them must then be assigned with a high RI score. The RT shift correction (§6.5)
/// on the same runs made 4 s late must put the reference Glycine 3TMS back on its library RI.
#[test]
fn amino_acids_get_ri_score_with_lib_aa() {
    use md_core::ri::{RiCalibration, marker_peaks, window_points};
    use md_core::search::{
        PreparedLibrary, SearchSettings, map_compounds_to_library, shift_to_reference,
    };
    let Ok(lbr) = std::env::var("MD_LBR_SAMPLE") else {
        eprintln!("skipped: MD_LBR_SAMPLE not set");
        return;
    };
    let runs = runs();
    let Some(dir) = runs.first().and_then(|p| p.parent()) else {
        return;
    };
    let lib = md_core::library::read_lbr(&lbr).unwrap();
    let lib_ri = |n: &str| lib.iter().find(|c| c.name == n).unwrap().spectrum.ri;
    let windows: Vec<(f64, f64, f64)> = table(&dir.join("rimLimits.txt"))
        .iter()
        .map(|r| {
            (
                r[0].parse().unwrap(),
                r[1].parse().unwrap(),
                r[2].parse().unwrap(),
            )
        })
        .collect();
    // (FAME RI in library.txt, name in lib_aa)
    let bridge = [(271500.0, "Valine 2TMS"), (368800.0, "Threonine 3TMS")];
    let mut best = [
        ("Leucine 2TMS", 0.0),
        ("Glycine 3TMS", 0.0),
        ("Serine 3TMS", 0.0),
    ];
    let s = SearchSettings {
        cutoff: 0.0,
        ..Default::default()
    };
    let prepared = PreparedLibrary::new(&lib, &s);
    let gly = lib.iter().position(|c| c.name == "Glycine 3TMS").unwrap();
    let mut corrected = 0;
    for path in &runs {
        let c = md_core::andi::read(path).unwrap();
        let fame = marker_peaks(&c, 87.0, DENSITY, THRESHOLD);
        let fame = RiCalibration::new(window_points(&fame, &windows)).unwrap();
        let cal = RiCalibration::new(bridge.map(|(f, n)| (fame.rt(f), lib_ri(n)))).unwrap();
        let mut traces = IonTraces::new(&c).unwrap();
        let peaks = traces.detect_peaks(&PeakSettings::default());
        let comps = deconvolute(&traces, &peaks, &c.scan_times_s, &DeconvSettings::default());
        let ris = cal.component_ris(&c.scan_times_s, &comps);
        let id = map_compounds_to_library(&comps, Some(&ris), &prepared, &s);
        for h in id.assigned.iter().flatten() {
            // Only clear spectrum matches count; weak ones are other substances.
            let (n, ri) = (&lib[h.compound].name, h.ri_sim.unwrap());
            if h.spec_sim > 0.65
                && let Some(b) = best.iter_mut().find(|b| b.0 == n)
            {
                eprintln!("{}: {n} spectrum {:.3}, RI {ri:.3}", name(path), h.spec_sim);
                b.1 = f64::max(b.1, ri);
            }
        }
        let late: Vec<f64> = c.scan_times_s.iter().map(|t| t + 4.0).collect();
        if let Some(shifted) = shift_to_reference(&cal, &late, &comps, &prepared, gly, 0.8, &s) {
            // The reference found is the component the search assigned to Glycine 3TMS.
            let k = id
                .assigned
                .iter()
                .position(|h| h.is_some_and(|h| h.compound == gly));
            let ri = shifted.ri(late[comps[k.unwrap()].apex_scan as usize]);
            eprintln!(
                "{}: shift {:.2} s, Glycine 3TMS at RI {ri:.2}",
                name(path),
                shifted.shift_s
            );
            assert!((ri - lib_ri("Glycine 3TMS")).abs() < 1e-6, "{}", name(path));
            let on_time =
                shift_to_reference(&cal, &c.scan_times_s, &comps, &prepared, gly, 0.8, &s);
            assert!((shifted.shift_s - on_time.unwrap().shift_s - 4.0).abs() < 1e-9);
            // Searched with the shifted RIs, the reference is assigned with its full RI score.
            let ris = shifted.component_ris(&late, &comps);
            let id = map_compounds_to_library(&comps, Some(&ris), &prepared, &s);
            let h = id.assigned[k.unwrap()].unwrap();
            assert_eq!(
                (h.compound, h.ri_sim.map(|r| r > 0.999_999)),
                (gly, Some(true))
            );
            corrected += 1;
        }
    }
    if !runs.is_empty() {
        // Each is found in at least one run, less than 0.46 sigma (9 RI) off its library RI (measured 0.94-0.99).
        assert!(best.iter().all(|b| b.1 > 0.9), "{best:?}");
        assert!(corrected > 0);
    }
}

/// MoNA GC-MS export (NIST MSP, `MD_MONA_MSP`; mixed licences, never in the repo) with the FAME RI
/// calibration of the example runs. The FiehnLib part of the file has `Retention_index:` on the same
/// scale as TargetSearch's `library.txt` (FAME RI in ms), so the markers calibrate the runs directly.
/// Prints, per run and over all runs, what the search finds with and without the RI score; a hit is
/// "wrong by RI" when the library RI is known and more than 3 sigma from the component. The targets of
/// `library.txt` count as found when a component inside the target's RI window is assigned a MoNA
/// substance of the same name.
#[test]
fn mona_msp_finds_the_targets_with_ri() {
    use md_core::ri::{RiCalibration, marker_peaks, window_points};
    use md_core::search::{PreparedLibrary, SearchSettings, map_compounds_to_library};
    let Ok(msp) = std::env::var("MD_MONA_MSP") else {
        eprintln!("skipped: MD_MONA_MSP not set");
        return;
    };
    let runs = runs();
    let Some(dir) = runs.first().and_then(|p| p.parent()) else {
        return;
    };
    let lib = md_core::library::import_msl(&msp).unwrap();
    let known = lib.iter().filter(|c| c.spectrum.ri >= 0.0).count();
    eprintln!(
        "MoNA: {} spectra, {known} with a retention index",
        lib.len()
    );
    let windows: Vec<(f64, f64, f64)> = table(&dir.join("rimLimits.txt"))
        .iter()
        .map(|r| {
            (
                r[0].parse().unwrap(),
                r[1].parse().unwrap(),
                r[2].parse().unwrap(),
            )
        })
        .collect();
    // (name without the derivative, RI, window)
    let targets: Vec<(String, f64, f64)> = table(&dir.join("library.txt"))
        .iter()
        .map(|r| {
            let name = r[0].split(" (").next().unwrap().to_lowercase();
            (name, r[1].parse().unwrap(), r[2].parse().unwrap())
        })
        .collect();
    // FAME RI is in ms: the sigma of the RI score is 2 s, not the 20 RI units of alkane indices.
    let with = SearchSettings {
        ri_diff: 2000.0,
        ..Default::default()
    };
    let without = SearchSettings {
        use_ri: false,
        ..with.clone()
    };
    let (prep_with, prep_without) = (
        PreparedLibrary::new(&lib, &with),
        PreparedLibrary::new(&lib, &without),
    );
    // [hits, wrong by RI, targets found] for [spectrum only, with RI]
    let mut total = [[0usize; 3]; 2];
    for path in &runs {
        let c = md_core::andi::read(path).unwrap();
        let fame = marker_peaks(&c, 87.0, DENSITY, THRESHOLD);
        let cal = RiCalibration::new(window_points(&fame, &windows)).unwrap();
        let mut traces = IonTraces::new(&c).unwrap();
        let peaks = traces.detect_peaks(&PeakSettings::default());
        let comps = deconvolute(&traces, &peaks, &c.scan_times_s, &DeconvSettings::default());
        let ris = cal.component_ris(&c.scan_times_s, &comps);
        let mut line = String::new();
        for (k, (s, prep)) in [(&without, &prep_without), (&with, &prep_with)]
            .into_iter()
            .enumerate()
        {
            let id = map_compounds_to_library(&comps, Some(&ris), prep, s);
            let mut n = [0usize; 3];
            let mut found = Vec::new();
            for (i, h) in id.assigned.iter().enumerate() {
                let Some(h) = h else { continue };
                let (c, ri) = (&lib[h.compound], ris[i]);
                n[0] += 1;
                n[1] +=
                    (c.spectrum.ri >= 0.0 && (c.spectrum.ri - ri).abs() > 3.0 * s.ri_diff) as usize;
                let name = c.name.to_lowercase();
                if let Some(t) = targets
                    .iter()
                    .find(|t| t.0 == name && (t.1 - ri).abs() <= t.2)
                {
                    n[2] += 1;
                    found.push(t.0.clone());
                }
            }
            line += &format!(
                " | {} {:?}: {found:?}",
                if k == 0 { "spec" } else { "RI" },
                &n[..2]
            );
            for (t, v) in total[k].iter_mut().zip(n) {
                *t += v;
            }
        }
        eprintln!("{}{line}", name(path));
    }
    eprintln!(
        "spectrum only: hits, wrong by RI, targets found = {:?}",
        total[0]
    );
    eprintln!(
        "with RI:       hits, wrong by RI, targets found = {:?}",
        total[1]
    );
    if !runs.is_empty() {
        // By construction the RI score drops every hit counted as wrong; the point is that the
        // spectrum alone does not.
        assert!(
            total[0][1] > total[1][1],
            "the RI score removed no wrong hit"
        );
        assert!(total[1][0] <= total[0][0]);
        assert!(total[1][2] >= total[0][2], "RI lost targets");
    }
}

/// #242: TargetSearch's own `library.txt` against the example runs (FAME calibration from the markers).
/// These runs have 20 scans/s: the ion apexes of one substance spread over several scans, so the default
/// deconvolution width of 1 scan cuts it into pieces below 10 ions and a run keeps ~15 components,
/// mostly the FAME markers (the old 3.3 gives the same, `old-results/`). A width of 10 scans and 5 ions
/// find the targets, and so does the automatic width of the defaults (#279: 10 scans here); the RI
/// sigma is on the library's FAME scale (2 s), not the alkane 20.
#[test]
fn targetsearch_library_finds_its_targets() {
    use md_core::ri::{RiCalibration, marker_peaks, window_points};
    use md_core::search::{PreparedLibrary, SearchSettings, map_compounds_to_library};
    let runs = runs();
    let Some(dir) = runs.first().and_then(|p| p.parent()) else {
        return;
    };
    let lib = md_core::library::import_targetsearch(dir.join("library.txt")).unwrap();
    let windows: Vec<(f64, f64, f64)> = table(&dir.join("rimLimits.txt"))
        .iter()
        .map(|r| {
            (
                r[0].parse().unwrap(),
                r[1].parse().unwrap(),
                r[2].parse().unwrap(),
            )
        })
        .collect();
    let s = SearchSettings {
        ri_diff: 2000.0,
        cutoff: 0.0,
        ..Default::default()
    };
    let prepared = PreparedLibrary::new(&lib, &s);
    let configs = [
        ("3.3 defaults (width 1 scan)", DeconvSettings::default()),
        (
            "defaults (automatic width, #279)",
            md_core::settings::GcmsSettings::default().deconv(),
        ),
        (
            "width 10, 5 ions",
            DeconvSettings {
                width: Some(10.0),
                min_peak_number: 5,
                ..Default::default()
            },
        ),
    ];
    // Runs in which each target is assigned with at least the batch cutoff 0.8, per config. 0.8 needs
    // an RI score above 0.51, i.e. within 2.3 s of the library RI.
    let mut found = configs.clone().map(|_| vec![0usize; lib.len()]);
    for path in &runs {
        let c = md_core::andi::read(path).unwrap();
        let cal = RiCalibration::new(window_points(
            &marker_peaks(&c, 87.0, DENSITY, THRESHOLD),
            &windows,
        ))
        .unwrap();
        let mut traces = IonTraces::new(&c).unwrap();
        let peaks = traces.detect_peaks(&PeakSettings::default());
        for (n, (_, d)) in configs.iter().enumerate() {
            let comps = deconvolute(&traces, &peaks, &c.scan_times_s, d);
            let ris = cal.component_ris(&c.scan_times_s, &comps);
            let id = map_compounds_to_library(&comps, Some(&ris), &prepared, &s);
            for (i, h) in id.assigned.iter().enumerate() {
                let Some(h) = h.filter(|h| h.overall >= 0.8) else {
                    continue;
                };
                found[n][h.compound] += 1;
                // A clear hit is also the first of the component's hit list.
                if n > 0 {
                    assert_eq!(id.hits[i][0].compound, h.compound, "{}", name(path));
                }
            }
        }
    }
    let targets = found
        .each_ref()
        .map(|f| f.iter().filter(|n| **n > 0).count());
    for ((what, _), (f, t)) in configs.iter().zip(found.iter().zip(targets)) {
        let per_target: Vec<_> = lib.iter().map(|c| c.name.as_str()).zip(f).collect();
        eprintln!("{what}: {t} targets, runs per target {per_target:?}");
    }
    let [old, auto, tuned] = targets;
    assert!(tuned >= 9, "{found:?}");
    assert!(old < tuned);
    // #279: the defaults find at least as many targets as the hand-tuned width.
    assert!(auto >= tuned, "automatic width {auto}, width 10 {tuned}");
}

/// Result table of a real run: components whose own spectra form the library are found again, with
/// their quant-ion area, and the CSV has one line per component.
#[test]
fn result_table_of_a_real_run() {
    use md_core::library::LibraryCompound;
    use md_core::results::{find_substance, result_rows, write_csv};
    use md_core::search::{PreparedLibrary, SearchSettings, map_compounds_to_library};

    let Some(path) = runs().into_iter().next() else {
        return;
    };
    let c = md_core::andi::read(&path).unwrap();
    let mut traces = IonTraces::new(&c).unwrap();
    let peaks = traces.detect_peaks(&PeakSettings::default());
    let masses: Vec<u32> = peaks
        .iter()
        .map(|p| traces.masses[p.row as usize])
        .collect();
    let comps = deconvolute(&traces, &peaks, &c.scan_times_s, &DeconvSettings::default());
    assert!(comps.len() > 10);
    let lib: Vec<LibraryCompound> = comps
        .iter()
        .take(3)
        .enumerate()
        .map(|(i, k)| {
            let mut l = LibraryCompound::new(&format!("Substance {i}"));
            l.spectrum.ri = -1.0;
            l.spectrum.masses = k.spectrum.iter().map(|s| s.0 as f32).collect();
            l.spectrum.values = k.spectrum.iter().map(|s| s.1 as i32).collect();
            l
        })
        .collect();
    let s = SearchSettings::default();
    let id = map_compounds_to_library(&comps, None, &PreparedLibrary::new(&lib, &s), &s);
    let rows = result_rows(&comps, &c.scan_times_s, None, &masses, &id, &lib, s.cutoff);
    assert_eq!(rows.len(), comps.len());
    for (i, r) in rows.iter().take(3).enumerate() {
        assert_eq!(r.name, format!("Substance {i}"));
        assert!(r.overall.unwrap() > 0.99 && r.area.unwrap() > 0.0, "{r:?}");
        assert_eq!(find_substance(&id, &lib, i, s.cutoff)[0].0, i);
    }
    let mut csv = Vec::new();
    write_csv(&mut csv, &rows).unwrap();
    assert_eq!(csv.iter().filter(|b| **b == b'\n').count(), rows.len() + 1);
}

/// MID wizard by sum formula (F16) on the example runs, which carry no label: the Glycine 3TMS
/// fragment m/z 174 (C7H20NSi2) of `lib_aa` must come out almost all M+0 once the natural isotopes
/// are taken out. Runs the whole pipeline: MID library, batch scan with rescan, `determine`.
#[test]
fn unlabelled_runs_give_m0_by_sum_formula() {
    use md_core::batch::{QuantSettings, ScanSettings, complete_scan, scan_runs};
    use md_core::isotope::IsotopeTable;
    use md_core::mid::{Method, determine, mid_library};
    use md_core::search::PreparedLibrary;
    use std::sync::atomic::{AtomicBool, AtomicUsize};
    let Ok(lbr) = std::env::var("MD_LBR_SAMPLE") else {
        eprintln!("skipped: MD_LBR_SAMPLE not set");
        return;
    };
    let runs = runs();
    if runs.is_empty() {
        return;
    }
    let mut lib = md_core::library::read_lbr(&lbr).unwrap();
    let gly = lib.iter().position(|c| c.name == "Glycine 3TMS").unwrap();
    lib[gly].formulas.push((174, 177, "C7H20NSi2".into()));
    let lib = mid_library(&lib);
    assert_eq!(lib.len(), 1);
    let mut s = ScanSettings::default();
    // Without RI calibration the spectra alone decide; the weaker matches still are Glycine.
    s.search.cutoff = 0.6;
    let prepared = PreparedLibrary::new(&lib, &s.search);
    let no = AtomicBool::new(false);
    let (mut found, errors) =
        scan_runs(&runs, &prepared, &s, 4, &AtomicUsize::new(0), &no).unwrap();
    assert!(errors.is_empty(), "{errors:?}");
    // Few runs find the substance without RI calibration; rescan the fragment ions in the others.
    let q = QuantSettings::default();
    complete_scan(&mut found, &lib, &s, &q, true, 4, &no).unwrap();
    let m = determine(
        &found,
        &lib,
        &[],
        &[],
        Method::Formula,
        "C",
        &IsotopeTable::standard(),
    );
    assert!(m.errors.is_empty(), "{:?}", m.errors);
    let fits: Vec<_> = m.compounds[0].fragments[0].fits.iter().flatten().collect();
    eprintln!(
        "Glycine 3TMS found with MIDs in {} of {} runs",
        fits.len(),
        runs.len()
    );
    assert_eq!(fits.len(), runs.len());
    for f in fits {
        eprintln!("M0..M3 {:.4?}, R² {:.4?}", f.mids, f.r2);
        assert!((f.mids[0] - 1.0).abs() < 0.03, "{:?}", f.mids);
        assert!(f.mids[1..].iter().all(|m| m.abs() < 0.03), "{:?}", f.mids);
    }
}

#[test]
fn component_detail_of_a_real_run() {
    use md_core::detail::{compare, component_detail};
    use md_core::library::LibraryCompound;
    use md_core::spectrum::Normalization;

    let Some(path) = runs().into_iter().next() else {
        return;
    };
    let c = md_core::andi::read(&path).unwrap();
    let mut traces = IonTraces::new(&c).unwrap();
    let peaks = traces.detect_peaks(&PeakSettings::default());
    let masses: Vec<u32> = peaks
        .iter()
        .map(|p| traces.masses[p.row as usize])
        .collect();
    let comps = deconvolute(&traces, &peaks, &c.scan_times_s, &DeconvSettings::default());
    let k = comps
        .iter()
        .find(|k| k.peaks.len() > 3)
        .expect("a component with ions");
    let d = component_detail(k, &peaks, &masses);
    assert!(d.peaks.len() <= k.peaks.len() && d.peaks.len() >= 3);
    assert!(d.mean_signal_to_noise > 0.0 && d.mean_sharpness > 0.0 && d.discrepancy >= 0.0);
    let mut lib = LibraryCompound::new("Self");
    lib.spectrum.masses = k.spectrum.iter().map(|s| s.0 as f32).collect();
    lib.spectrum.values = k.spectrum.iter().map(|s| s.1 as i32).collect();
    let cmp = compare(&k.spectrum, &lib.spectrum, Normalization::BasePeak);
    assert_eq!(cmp.found, cmp.library);
    assert!(cmp.difference().iter().all(|d| *d == 0.0));
}

/// Saves the components of two example runs as `.cmp` into `example-data/old-results/md-core/` (ignored
/// by git; `make-reference.sh --reread` checks there that the old program reads them) and loads them
/// again with the values the component panel shows.
#[test]
fn example_components_save_and_load_again() {
    use md_core::cmp::{components_from_file, compounds_of_components, write_cmp};
    use md_core::detail::component_detail;

    let out = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../example-data/old-results/md-core");
    for path in runs().into_iter().take(2) {
        let c = md_core::andi::read(&path).unwrap();
        let mut traces = IonTraces::new(&c).unwrap();
        let peaks = traces.detect_peaks(&PeakSettings::default());
        let masses: Vec<u32> = peaks
            .iter()
            .map(|p| traces.masses[p.row as usize])
            .collect();
        let mut comps = deconvolute(&traces, &peaks, &c.scan_times_s, &DeconvSettings::default());
        comps.sort_by_key(|k| k.apex_scan);
        assert!(comps.len() > 10, "{}", name(&path));

        let saved = compounds_of_components(&comps, &masses, &peaks, &c.scan_times_s, &[]);
        std::fs::create_dir_all(&out).unwrap();
        let file = out.join(Path::new(&*name(&path)).with_extension("cmp"));
        write_cmp(&file, &saved).unwrap();
        let (back, m2, p2) = components_from_file(&file, &c.scan_times_s).unwrap();

        assert_eq!(back.len(), comps.len());
        for (a, b) in comps.iter().zip(&back) {
            assert_eq!(
                (
                    a.first_scan,
                    a.apex_scan,
                    a.last_scan,
                    &a.model,
                    &a.spectrum
                ),
                (
                    b.first_scan,
                    b.apex_scan,
                    b.last_scan,
                    &b.model,
                    &b.spectrum
                )
            );
            // The file keeps one peak per mass in ascending order; the panel values are the same.
            let (mut before, after) = (
                component_detail(a, &peaks, &masses),
                component_detail(b, &p2, &m2),
            );
            before.peaks.sort_by_key(|r| r.mass);
            assert_eq!(before, after, "{} at scan {}", name(&path), a.apex_scan);
        }
    }
}
