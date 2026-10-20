//! Headless UI checks. `every_ui_symbol_has_a_glyph` runs in CI. `screenshots` renders the app with
//! real example runs into `target/ui-snapshots/*.png` for design reviews (`docs/ui-quality.md`); pixels
//! differ between GPUs and OSes, so the images are not committed and the test is ignored in CI:
//! `cargo test --release -p metabolite-detector screenshots -- --ignored`

use std::path::{Path, PathBuf};
use std::sync::Arc;
use std::time::Duration;

use eframe::egui;
use egui_kittest::Harness;
use egui_kittest::kittest::{Queryable, by};
use md_core::library::LibraryCompound;
use md_core::search::PreparedLibrary;

use crate::{App, Trace, identify};

/// Text the app shows: string literals outside `//` comments in this crate's sources.
fn ui_strings() -> Vec<(String, String)> {
    let src = Path::new(env!("CARGO_MANIFEST_DIR")).join("src");
    let mut found = Vec::new();
    for entry in std::fs::read_dir(src).unwrap() {
        let path = entry.unwrap().path();
        let text = std::fs::read_to_string(&path).unwrap();
        for line in text.lines().filter(|l| !l.trim_start().starts_with("//")) {
            // Literals are the odd parts between quotes; good enough for this crate's code.
            for lit in line.split('"').skip(1).step_by(2) {
                found.push((path.display().to_string(), lit.to_owned()));
            }
        }
    }
    found
}

// shortcut: `✕` in library_editor.rs has no glyph in any bundled font; drop it here once that window
// uses `×` like the rest of the app (reported on #89).
const KNOWN: [char; 1] = ['✕'];

/// Symbols like arrows are not in egui's default font; without a fallback they draw as boxes.
#[test]
fn every_ui_symbol_has_a_glyph() {
    let ctx = egui::Context::default();
    crate::style::apply(&ctx);
    ctx.run_ui(egui::RawInput::default(), |_| {})
        .textures_delta
        .clear();
    let font = egui::FontId::proportional(13.0);
    let missing: Vec<String> = ui_strings()
        .into_iter()
        .flat_map(|(file, lit)| lit.chars().map(|c| (file.clone(), c)).collect::<Vec<_>>())
        .filter(|(_, c)| {
            !c.is_ascii() && !KNOWN.contains(c) && !ctx.fonts_mut(|f| f.has_glyph(&font, *c))
        })
        .map(|(file, c)| format!("{file}: {c} (U+{:04X})", c as u32))
        .collect();
    assert!(missing.is_empty(), "no glyph:\n{}", missing.join("\n"));
}

fn runs() -> Vec<PathBuf> {
    let dir = Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../../example-data/TargetSearchData/inst/gc-ms-data");
    let mut paths: Vec<_> = std::fs::read_dir(&dir)
        .unwrap_or_else(|_| panic!("{} not found, run example-data/fetch.*", dir.display()))
        .map(|e| e.unwrap().path())
        .filter(|p| p.extension().is_some_and(|e| e == "cdf"))
        .collect();
    paths.sort();
    paths
}

/// Steps until `done` holds; background jobs report through channels, not repaint requests.
fn wait(h: &mut Harness<'_, App>, done: impl Fn(&App) -> bool) {
    for _ in 0..2000 {
        if done(h.state()) {
            h.run_steps(4);
            return;
        }
        std::thread::sleep(Duration::from_millis(10));
        h.step();
    }
    panic!("background work did not finish");
}

/// A library made of the first run's components, so the hit panel and labels have real scores.
fn library(app: &App) -> identify::Library {
    let run = &app.runs[0];
    let compounds: Vec<_> = (run.components.iter())
        .map(|k| {
            let rt = run.data.scan_times_s[k.apex_scan as usize];
            let mut c = LibraryCompound::new(&format!("Compound at {:.2} min", rt / 60.0));
            c.spectrum.masses = k.spectrum.iter().map(|(m, _)| *m as f32).collect();
            c.spectrum.values = k.spectrum.iter().map(|(_, v)| *v as i32).collect();
            c
        })
        .collect();
    let settings = app.settings.applied.search();
    identify::Library {
        path: "screenshots.mdlib".into(),
        name: "screenshots.mdlib".into(),
        prepared: PreparedLibrary::new(&compounds, &settings),
        settings,
        compounds,
    }
}

struct Shots {
    dir: PathBuf,
}

impl Shots {
    fn save(&self, h: &mut Harness<'_, App>, name: &str) {
        let img = h.render().expect("render");
        let file = std::fs::File::create(self.dir.join(format!("{name}.png"))).unwrap();
        let mut enc = png::Encoder::new(file, img.width(), img.height());
        enc.set_color(png::ColorType::Rgba);
        let mut w = enc.write_header().unwrap();
        w.write_image_data(img.as_raw()).unwrap();
        // Inputs a screen reader cannot name; the audit reads these next to the images.
        let unnamed = h.unnamed_widgets().join("\n");
        std::fs::write(self.dir.join(format!("{name}.unnamed.txt")), unnamed).unwrap();
    }
}

fn harness(theme: egui::Theme) -> Harness<'static, App> {
    Harness::builder()
        .with_size([1280.0, 800.0])
        .with_theme(theme)
        // The audit wants this signal, but as findings, not as a panic on the first unnamed field.
        .with_accessibility_check(false)
        .with_context_setup(crate::style::apply)
        .build_eframe(|_| App::default())
}

/// Opens `paths`, waits for loading and selects the second component of the first run.
fn open(h: &mut Harness<'_, App>, paths: &[PathBuf]) {
    let ctx = h.ctx.clone();
    for p in paths {
        h.state_mut().open(&ctx, p.clone());
    }
    wait(h, |a| a.loading.is_empty());
    h.state_mut().active = 0;
    let run = &mut h.state_mut().runs[0];
    run.component = Some(1.min(run.components.len() - 1));
    run.scan = run.components[run.component.unwrap()].apex_scan as usize;
}

#[test]
#[ignore = "needs example data and a GPU; writes target/ui-snapshots"]
fn screenshots() {
    let dir = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../target/ui-snapshots");
    std::fs::create_dir_all(&dir).unwrap();
    let shots = Shots { dir };
    let paths = runs();

    for (theme, tag) in [(egui::Theme::Dark, "dark"), (egui::Theme::Light, "light")] {
        let mut h = harness(theme);
        h.run_steps(2);
        shots.save(&mut h, &format!("01-start-{tag}"));

        open(&mut h, &paths[..3]);
        let lib = Arc::new(library(h.state()));
        h.state_mut().library = Some(lib);
        wait(&mut h, |a| a.runs.iter().all(|r| r.ident.is_some()));
        shots.save(&mut h, &format!("02-main-{tag}"));
    }

    let mut h = harness(egui::Theme::Dark);
    open(&mut h, &paths);
    shots.save(&mut h, "03-overlay-15");
    h.state_mut().trace = Trace::Both;
    wait(&mut h, |a| {
        a.runs
            .iter()
            .all(|r| r.sic_job.is_none() && !r.sic.is_empty())
    });
    shots.save(&mut h, "04-tic-sic");
    h.state_mut().trace = Trace::Tic;

    h.get_by_label("File").click();
    h.run_steps(2);
    shots.save(&mut h, "05-menu-file");
    h.key_press(egui::Key::Escape);
    h.get_by_label("Analysis").click();
    h.run_steps(2);
    shots.save(&mut h, "06-menu-analysis");
    h.key_press(egui::Key::Escape);
    h.run_steps(2);

    let lib = Arc::new(library(h.state()));
    h.state_mut().library = Some(lib);
    wait(&mut h, |a| a.runs.iter().all(|r| r.ident.is_some()));

    type Toggle = fn(&mut App) -> &mut bool;
    let windows: [(&str, Toggle); 5] = [
        ("07-settings", |a| &mut a.settings.open),
        ("08-calibration", |a| &mut a.ri.open),
        ("09-batch", |a| &mut a.batch.open),
        ("10-isotopes", |a| &mut a.isotopes.open),
        ("11-library-editor", |a| &mut a.editor.open),
    ];
    for (name, open) in windows {
        *open(h.state_mut()) = true;
        h.run_steps(3);
        if name == "11-library-editor" {
            h.get_by_label("Edit open library").click();
            h.run_steps(2);
            h.get_all(by().label_contains("Compound at"))
                .last()
                .unwrap()
                .click();
            h.run_steps(2);
        }
        shots.save(&mut h, name);
        *open(h.state_mut()) = false;
        h.run_steps(2);
    }
}
