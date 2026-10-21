//! Headless UI checks. `every_ui_symbol_has_a_glyph` runs in CI. `screenshots` renders the app with
//! real example runs into `target/ui-snapshots/*.png` (each also as `*.de.png` in German) for design reviews (`docs/ui-quality.md`); pixels
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

use crate::guided::{Flow, Mode};
use crate::{App, Trace, identify};

/// Non-ASCII characters on code lines (not `//` comments) of this crate's sources. Code outside
/// strings and char literals is ASCII, so these are the symbols the app can show.
fn ui_symbols() -> Vec<(String, char)> {
    let mut found = Vec::new();
    for (path, text) in sources() {
        for line in text.lines().filter(|l| !l.trim_start().starts_with("//")) {
            let chars = line.chars().filter(|c| !c.is_ascii());
            found.extend(chars.map(|c| (path.display().to_string(), c)));
        }
    }
    found
}

/// This crate's `.rs` files and their text.
fn sources() -> Vec<(PathBuf, String)> {
    let mut dirs = vec![Path::new(env!("CARGO_MANIFEST_DIR")).join("src")];
    let mut found = Vec::new();
    while let Some(dir) = dirs.pop() {
        for entry in std::fs::read_dir(dir).unwrap() {
            let path = entry.unwrap().path();
            if path.is_dir() {
                dirs.push(path);
            } else if path.extension().is_some_and(|e| e == "rs") {
                let text = std::fs::read_to_string(&path).unwrap();
                found.push((path, text));
            }
        }
    }
    found
}

/// Calls whose first string argument is shown to the user.
const TEXT_CALLS: [&str; 27] = [
    "label",
    "button",
    "heading",
    "strong",
    "small",
    "weak",
    "on_hover_text",
    "on_disabled_hover_text",
    "hint_text",
    "menu_button",
    "selectable_label",
    "selectable_value",
    "checkbox",
    "radio",
    "radio_value",
    "toggle_value",
    "hyperlink_to",
    "text",
    "add_filter",
    "set_title",
    "set_description",
    "Window::new",
    "CollapsingHeader::new",
    "RichText::new",
    "Button::new",
    "Label::new",
    "show_note",
];

/// Whether a literal is prose to translate: it has a lowercase word ("Open", "Läufe"), not only
/// symbols, units or abbreviations ("RI", "m/z", "S/N"). Format placeholders do not count.
fn is_prose(literal: &str) -> bool {
    let mut plain = String::new();
    let mut depth = 0;
    for c in literal.chars() {
        match c {
            '{' => depth += 1,
            '}' => depth -= 1,
            _ if depth == 0 => plain.push(c),
            _ => {}
        }
    }
    // A word: two letters or more, all lowercase after the first ("Run", "ms"; not "RI", "NetCDF").
    (plain.split(|c: char| !c.is_alphabetic())).any(|w| {
        let mut rest = w.chars().skip(1).peekable();
        rest.peek().is_some() && rest.all(char::is_lowercase)
    })
}

/// String literals of the argument list `args` (the text after a call's `(`) that are passed as they
/// are: directly, or as the format string of a `format!`. Literals inside other nested calls, like
/// `t(..)`, or inside `[..]` (file extensions) are not.
/// shortcut: char literals like `'('` are not skipped; none sits in a UI call today.
fn bare_literals(args: &str) -> Vec<&str> {
    let (mut found, mut depth, mut format_at) = (Vec::new(), 0, None);
    let mut chars = args.char_indices();
    while let Some((i, c)) = chars.next() {
        match c {
            '"' => {
                let mut end = args.len();
                let mut escaped = false;
                for (j, c) in chars.by_ref() {
                    match c {
                        '"' if !escaped => {
                            end = j;
                            break;
                        }
                        '\\' => escaped = !escaped,
                        _ => escaped = false,
                    }
                }
                if depth == 0 || format_at == Some(depth) {
                    found.push(&args[i + 1..end]);
                }
            }
            '(' | '[' => {
                depth += 1;
                if depth == 1 && c == '(' && args[..i].trim_end().ends_with("format!") {
                    format_at = Some(1);
                }
            }
            ')' | ']' if depth == 0 => break,
            ')' | ']' => {
                if format_at == Some(depth) {
                    format_at = None;
                }
                depth -= 1;
            }
            _ => {}
        }
    }
    found
}

/// Texts passed to a [`TEXT_CALLS`] call as a bare literal or `format!`, not through `lang::t` or
/// `lang::tf!`. Test modules are skipped.
/// shortcut: texts that reach the UI another way (a const array, a variable) are not seen; the screenshot
/// check in German catches those.
fn untranslated() -> Vec<String> {
    let mut found = Vec::new();
    for (path, text) in sources() {
        let name = path.file_name().unwrap().to_string_lossy();
        if name == "ui_tests.rs" || name == "lang.rs" {
            continue;
        }
        let text = text.replace("\r\n", "\n");
        let text = &text[..text.find("#[cfg(test)]\nmod ").unwrap_or(text.len())];
        let code: String = (text.lines())
            .map(|l| {
                if l.trim_start().starts_with("//") {
                    ""
                } else {
                    l
                }
            })
            .collect::<Vec<_>>()
            .join("\n");
        for call in TEXT_CALLS {
            let pattern = format!("{call}(");
            for (at, _) in code.match_indices(&pattern) {
                let before = code[..at].chars().next_back();
                // `plot_export::button(..)` is our own function, not the widget method.
                let path_call = before == Some(':') && !call.contains("::");
                if path_call || before.is_some_and(|c| c.is_alphanumeric() || c == '_') {
                    continue;
                }
                // The first string is the shown text; later ones are file names or keys.
                for literal in bare_literals(&code[at + pattern.len()..])
                    .into_iter()
                    .take(1)
                {
                    if is_prose(literal) {
                        let line = code[..at].lines().count();
                        found.push(format!("{}:{line}: {call}(\"{literal}\"", path.display()));
                    }
                }
            }
        }
    }
    found
}

/// S8: every UI text has an English and a German version (#257).
#[test]
fn every_ui_text_is_translated() {
    let missing = untranslated();
    assert!(
        missing.is_empty(),
        "{} texts without lang::t or lang::tf!:\n{}",
        missing.len(),
        missing.join("\n")
    );
}

/// S8: an error shown in red is a fixed text or goes through `errors::text`, which translates what
/// `md-core` says (#339). `text` is the hit list of the batch table, not an error.
#[test]
fn error_labels_translate_core_messages() {
    let mut found = Vec::new();
    for (path, text) in sources() {
        if path.file_name().is_some_and(|n| n == "ui_tests.rs") {
            continue;
        }
        let text = text.replace("\r\n", "\n");
        let text = &text[..text.find("#[cfg(test)]\nmod ").unwrap_or(text.len())];
        for key in [
            "colored_label(ui.visuals().error_fg_color,",
            "colored_label(error,",
        ] {
            for (at, _) in text.match_indices(key) {
                let arg = text[at + key.len()..].trim_start();
                let ok = ["errors::text(", "t(", "tf!(", "lang::t(", "text)"]
                    .iter()
                    .any(|p| arg.starts_with(p));
                if !ok {
                    let line = text[..at].lines().count() + 1;
                    found.push(format!("{}:{line}", path.display()));
                }
            }
        }
    }
    assert!(
        found.is_empty(),
        "error labels without errors::text:\n{}",
        found.join("\n")
    );
}

#[test]
fn prose_needs_a_lowercase_word() {
    assert!(is_prose("Open a run"));
    assert!(is_prose("{n} runs"));
    assert!(!is_prose("RI"));
    assert!(!is_prose("m/z {mz:.1}"));
    assert!(is_prose("Run") && is_prose("New"));
    assert!(!is_prose("NetCDF (ANDI-MS)"));
}

#[test]
fn bare_literals_skip_translated_and_nested_ones() {
    let args = |s| bare_literals(s);
    assert_eq!(
        args(r#"x.is_some(), "Open results") + 1"#),
        ["Open results"]
    );
    assert_eq!(
        args(r#"t("Open", "Öffnen")).clicked()"#),
        Vec::<&str>::new()
    );
    assert_eq!(args(r#"format!("{n} runs", n = f(1)))"#), ["{n} runs"]);
    assert_eq!(args(r#"t("CSV", "CSV"), &["csv"])"#), Vec::<&str>::new());
    assert_eq!(args(r#""say \"hi\"")"#), [r#"say \"hi\""#]);
}

/// Symbols like arrows are not in egui's default font; without a fallback they draw as boxes.
#[test]
fn every_ui_symbol_has_a_glyph() {
    let ctx = egui::Context::default();
    crate::style::apply(&ctx);
    ctx.run_ui(egui::RawInput::default(), |_| {})
        .textures_delta
        .clear();
    let font = egui::FontId::proportional(13.0);
    let missing: Vec<String> = ui_symbols()
        .into_iter()
        .filter(|(_, c)| !ctx.fonts_mut(|f| f.has_glyph(&font, *c)))
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
    let mut compounds: Vec<_> = (run.components.iter())
        .map(|k| {
            let rt = run.data.scan_times_s[k.apex_scan as usize];
            let mut c = LibraryCompound::new(&format!("Compound at {:.2} min", rt / 60.0));
            c.spectrum.masses = k.spectrum.iter().map(|(m, _)| *m as f32).collect();
            c.spectrum.values = k.spectrum.iter().map(|(_, v)| *v as i32).collect();
            c
        })
        .collect();
    // The first substance shows the fields of the newer library versions.
    if let Some(c) = compounds.first_mut() {
        c.auto_quant_ions = c.spectrum.masses.iter().take(3).copied().collect();
        c.ion_freq = c.spectrum.masses.iter().map(|&m| (m, 3)).collect();
        c.spec_count = 3;
        c.formulas.push((73, 76, "C3H9Si".into()));
        c.mid_sizes.push((73, 76, 4));
    }
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
    /// Saves `name.png` and, drawn with the same state in German, `name.de.png` (S8, #257). The test
    /// steers by English labels, so it switches only for the picture.
    fn save(&self, h: &mut Harness<'_, App>, name: &str) {
        self.write(h, name);
        crate::lang::set(crate::lang::Lang::De);
        h.run_steps(2);
        self.write(h, &format!("{name}.de"));
        crate::lang::set(crate::lang::Lang::En);
        h.run_steps(2);
    }

    fn write(&self, h: &mut Harness<'_, App>, name: &str) {
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
    // Images of a renamed or dropped screen would otherwise sit next to the new ones.
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();
    let shots = Shots { dir };
    let paths = runs();

    for (theme, tag) in [(egui::Theme::Dark, "dark"), (egui::Theme::Light, "light")] {
        let mut h = harness(theme);
        h.run_steps(2);
        shots.save(&mut h, &format!("01-start-{tag}"));
        // As `persist::restore` asks: Guided preselected.
        (h.state_mut().guided.mode, h.state_mut().guided.ask) = (Mode::Guided, true);
        h.run_steps(3);
        shots.save(&mut h, &format!("01b-first-start-{tag}"));
        (h.state_mut().guided.mode, h.state_mut().guided.ask) = (Mode::Expert, false);

        open(&mut h, &paths[..3]);
        shots.save(&mut h, &format!("02a-main-no-library-{tag}"));
        let lib = Arc::new(library(h.state()));
        h.state_mut().library = Some(lib);
        wait(&mut h, |a| a.runs.iter().all(|r| r.ident.is_some()));
        shots.save(&mut h, &format!("02-main-{tag}"));
        guided_shots(&shots, &mut h, tag);
        // The triangles on the time axis stay with the derivative and in the SIC.
        h.state_mut().view.derivative = true;
        h.state_mut().trace = Trace::Both;
        wait(&mut h, |a| {
            a.runs
                .iter()
                .all(|r| r.sic_job.is_none() && !r.sic.is_empty())
        });
        shots.save(&mut h, &format!("02c-derivative-sic-{tag}"));
    }
    // The smallest window the layout is checked at: the spectrum must stay in view.
    let mut h = Harness::builder()
        .with_size([800.0, 600.0])
        .with_accessibility_check(false)
        .with_context_setup(crate::style::apply)
        .build_eframe(|_| App::default());
    open(&mut h, &paths[..3]);
    h.state_mut().model_profile = true;
    h.run_steps(4);
    shots.save(&mut h, "02b-main-800x600");

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
    h.get_by_label("Export chromatogram").click();
    h.run_steps(2);
    shots.save(&mut h, "05b-menu-export-plot");
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
    let windows: [(&str, Toggle); 6] = [
        ("07-settings", |a| &mut a.settings.open),
        ("08-calibration", |a| &mut a.ri.open),
        ("09-batch", |a| &mut a.batch.open),
        ("10-isotopes", |a| &mut a.isotopes.open),
        ("11-library-editor", |a| &mut a.editor.open),
        ("13-run-info", |a| &mut a.view.info_open),
    ];
    for (name, open) in windows {
        *open(h.state_mut()) = true;
        h.run_steps(3);
        if name == "11-library-editor" {
            h.get_by_label("Edit open library").click();
            h.run_steps(2);
            // The run list and the editor's list both show "Compound at"; click the first entry of the
            // editor's list (its left column), not whatever the tree ends with.
            let win = h.get_by_label_contains("Library editor").rect();
            h.get_all(by().label_contains("Compound at"))
                .find(|n| win.contains_rect(n.rect()) && n.rect().min.x < win.min.x + 40.0)
                .expect("no compound in the library editor list")
                .click();
            h.run_steps(2);
            for section in ["Automatic quant ions", "MID fragments (sum formulas)"] {
                h.get_by_label(section).click();
                h.run_steps(2);
            }
            // A grid takes its column widths from the frame before.
            h.run_steps(3);
        }
        shots.save(&mut h, name);
        *open(h.state_mut()) = false;
        h.run_steps(2);
    }

    // MIDs of one substance with a made-up fragment: its four largest-mass neighbours from the top ion.
    let mut lib = library(h.state());
    let c = &mut lib.compounds[1];
    let top = (c.spectrum.masses.iter().zip(&c.spectrum.values))
        .max_by_key(|p| *p.1)
        .map_or(73, |p| *p.0 as u32);
    c.formulas.push((top, top + 3, "C7H20NSi2".into()));
    h.state_mut().library = Some(Arc::new(lib));
    h.state_mut().mid.demo(paths[..4].to_vec());
    h.state_mut().mid.open = true;
    h.run_steps(3);
    shots.save(&mut h, "14-mid-settings");
    h.get_by_label("Start").click();
    h.run_steps(2);
    wait(&mut h, |a| !a.mid.busy());
    h.get_all(by().label_contains("Compound at"))
        .next()
        .unwrap()
        .click();
    h.run_steps(3);
    shots.save(&mut h, "15-mid-results");
    h.state_mut().mid.open = false;

    // A batch of made-up quantities over the open runs, in groups of ten by file name, so the
    // table and the statistics have something to show without a real batch run.
    let lib = h.state().library.clone().unwrap();
    // The library is the first run's components, so compound k sits where component k does.
    let run = &h.state().runs[0];
    let rts: Vec<f64> = (run.components.iter())
        .map(|k| run.data.scan_times_s[k.apex_scan as usize])
        .collect();
    let found: Vec<_> = (paths.iter().enumerate())
        .map(|(i, p)| md_core::batch::RunResult {
            path: p.clone(),
            found: (0..lib.compounds.len().min(16))
                .map(|k| md_core::batch::Found {
                    compound: k,
                    score: 0.9,
                    rt_s: rts[k],
                    ions: [(73, 1.0), (147, 0.4)]
                        .map(|(mass, share)| md_core::batch::IonPeak {
                            mass,
                            integral: share * fake_integral(p, i, k),
                            quality: 0.1,
                            signal_to_noise: 10.0,
                        })
                        .into(),
                    rescanned: false,
                    rescued: Vec::new(),
                    alternatives: Vec::new(),
                })
                .collect(),
            ..Default::default()
        })
        .collect();
    let batch = &mut h.state_mut().batch;
    batch.finished(lib, found);
    batch.settings.replicates = r"eg(\d)\d".into();
    batch.open = true;
    h.run_steps(3);
    shots.save(&mut h, "16-batch-table");
    h.state_mut().batch.open = false;
    h.state_mut().statistics.open = true;
    h.run_steps(3);
    // Group 0 (runs 04 to 09) against group 2 (runs 20 to 26).
    h.state_mut().statistics.compare("0", "2");
    h.run_steps(3);
    shots.save(&mut h, "17-statistics");
    h.hover_at(egui::pos2(800.0, 500.0));
    h.run_steps(2);
    h.event(egui::Event::MouseWheel {
        unit: egui::MouseWheelUnit::Point,
        delta: egui::vec2(-600.0, 0.0),
        phase: egui::TouchPhase::Move,
        source: egui::MouseWheelSource::Wheel,
        modifiers: egui::Modifiers::NONE,
    });
    h.run_steps(5);
    shots.save(&mut h, "17b-statistics-scrolled");
}

/// Guided mode: the Settings window with "Advanced" closed and one step of each flow.
fn guided_shots(shots: &Shots, h: &mut Harness<'_, App>, tag: &str) {
    h.state_mut().guided.mode = Mode::Guided;
    h.state_mut().settings.open = true;
    h.run_steps(3);
    shots.save(h, &format!("18-guided-settings-{tag}"));
    h.state_mut().settings.open = false;
    for (flow, step, name) in [
        (Flow::Ri, 1, "19-guided-ri"),
        (Flow::Detection, 2, "20-guided-detection"),
        (Flow::Batch, 2, "21-guided-batch"),
        (Flow::Mid, 1, "22-guided-mid"),
    ] {
        h.state_mut().start_flow(flow);
        h.state_mut().guided.flow = Some((flow, step));
        h.run_steps(3);
        shots.save(h, &format!("{name}-{tag}"));
    }
    h.state_mut().guided.flow = None;
    h.state_mut().guided.mode = Mode::Expert;
    h.run_steps(2);
}

/// A quantity that differs by substance, a little by run and, for most substances, by the run's
/// group (the tens digit of its file name, as the replicate regex above takes it), so p-values
/// differ between rows.
fn fake_integral(path: &Path, run: usize, compound: usize) -> f64 {
    let name = path.file_stem().unwrap().to_string_lossy();
    let group = f64::from(name.as_bytes()[name.len() - 2] - b'0');
    let noise = ((run * 7 + compound * 3) % 5) as f64;
    let effect = 0.1 * (compound % 5) as f64;
    1000.0 * (compound + 1) as f64 * (1.0 + effect * group + 0.05 * noise)
}

#[test]
#[ignore = "needs a GPU; writes target/ui-snapshots"]
fn tof_calibration_screenshot() {
    let dir = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../target/ui-snapshots");
    std::fs::create_dir_all(&dir).unwrap();
    let shots = Shots { dir };
    for (theme, tag) in [(egui::Theme::Dark, "dark"), (egui::Theme::Light, "light")] {
        let mut h = harness(theme);
        h.state_mut().tof = crate::tof_calibration::TofWindow::demo();
        h.run_steps(3);
        shots.save(&mut h, &format!("14-tof-calibration-{tag}"));
    }
}

/// The view toolbar: a score triangle switches its band, the zoom buttons run, "?" explains the plot.
#[test]
fn the_view_toolbar_works_on_a_run() {
    let mut h = harness(egui::Theme::Dark);
    h.state_mut().runs.push(crate::results::tests::run());
    h.run_steps(3);
    h.get_by_label("▲ ≤ 0.7").click();
    h.get_all(by().label("+")).next().unwrap().click();
    h.run_steps(2);
    h.get_all(by().label("All")).next().unwrap().click();
    h.run_steps(2);
    assert_eq!(h.state().view.score_filter, [false, true, true]);
    h.get_all(by().label("?")).next().unwrap().click();
    h.run_steps(2);
    h.get_by_label_contains("Double-click, or All");
}

/// First start: Guided is preselected and the answer sticks; Guided closes "Advanced", Expert opens it;
/// a flow can be skipped through and ends in its expert window.
#[test]
fn guided_mode_asks_folds_advanced_and_ends_in_the_expert_window() {
    let mut h = harness(egui::Theme::Dark);
    (h.state_mut().guided.mode, h.state_mut().guided.ask) = (Mode::Guided, true);
    h.run_steps(2);
    h.get_by_label("Continue").click();
    h.run_steps(2);
    assert!(!h.state().guided.ask);
    assert_eq!(h.state().guided.mode, Mode::Guided);

    h.state_mut().settings.open = true;
    h.run_steps(2);
    h.get_by_label("Advanced");
    assert!(h.query_by_label("Peak detection method").is_none());
    h.state_mut().guided.mode = Mode::Expert;
    h.run_steps(2);
    h.get_by_label("Peak detection method");
    h.state_mut().settings.open = false;

    h.state_mut().start_flow(Flow::Ri);
    h.run_steps(2);
    h.get_by_label_contains("Step 1 of 3");
    for _ in 0..2 {
        h.get_by_label("Skip").click();
        h.run_steps(2);
    }
    h.get_by_label_contains("Step 3 of 3");
    h.get_by_label("Finish").click();
    h.run_steps(2);
    assert!(h.state().guided.flow.is_none());
    assert!(h.state().ri.open);

    // A detection switches .cmp saving on; Finish, cancel and another flow switch it off again.
    let app = h.state_mut();
    app.start_flow(Flow::Detection);
    assert!(app.batch.saves_cmp());
    app.start_flow(Flow::Batch);
    assert!(!app.batch.saves_cmp());
    app.start_flow(Flow::Detection);
    app.cancel_flow();
    assert!(!app.batch.saves_cmp());
    app.start_flow(Flow::Detection);
    app.guided.flow = Some((Flow::Detection, 3));
    h.run_steps(2);
    h.get_by_label("Finish").click();
    h.run_steps(2);
    assert!(h.state().batch.open && !h.state().batch.saves_cmp());
}
