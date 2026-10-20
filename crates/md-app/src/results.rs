//! Result table of the active run (old `DataTableModel` with `CompoundExporter` and Find Compound): one
//! row per component, sortable, a click selects the component, a substance can be looked up in the
//! library and the table exported as CSV.

use eframe::egui;
use egui_extras::{Column, TableBuilder};
use md_core::results::{ResultRow, find_substance, result_rows, write_csv};

use crate::Run;
use crate::identify::Library;

const HEAD: [&str; 8] = [
    "Compound Name",
    "RT (min)",
    "RI",
    "Peak Area",
    "Overall Score",
    "Spec. Score",
    "RI Score",
    "Quant. Ions",
];

/// What the cached rows were built from: active run and its `rev`.
type Stamp = (usize, u64);

#[derive(Default)]
pub struct Results {
    pub open: bool,
    /// Only rows whose name contains this (case-insensitive), and the library names offered for a lookup.
    filter: String,
    /// Column and descending.
    sort: Option<(usize, bool)>,
    /// Library compound whose components are listed.
    looked_up: Option<usize>,
    note: String,
    cache: Option<(Stamp, Vec<ResultRow>)>,
}

impl Results {
    pub fn window(&mut self, ctx: &egui::Context, runs: &mut [Run], active: usize) {
        let mut open = self.open;
        egui::Window::new("Results")
            .open(&mut open)
            .default_size([900.0, 500.0])
            .show(ctx, |ui| self.contents(ui, runs.get_mut(active), active));
        self.open = open;
    }

    fn contents(&mut self, ui: &mut egui::Ui, run: Option<&mut Run>, active: usize) {
        let Some(run) = run else {
            ui.label("Open a chromatogram first.");
            return;
        };
        let Some(((lib, ..), _, _)) = &run.ident else {
            ui.horizontal(|ui| self.saved_buttons(ui, run));
            if !self.note.is_empty() {
                ui.label(&self.note);
            }
            ui.label("Open a library to identify the components of the run.");
            return;
        };
        let lib = lib.clone();
        let stamp = (active, run.rev);
        let mut cache = self.cache.take().filter(|c| c.0 == stamp);
        if cache.is_none() {
            cache = Some((stamp, build_rows(run)));
            self.looked_up = None;
        }
        let (stamp, rows) = cache.expect("filled above");

        let mut loaded = false;
        ui.horizontal(|ui| {
            ui.label(format!("{} — library {}", run.name, lib.name));
            if ui.button("Export CSV…").clicked() {
                self.note = export(&run.name, &rows);
            }
            loaded = self.saved_buttons(ui, run);
        });
        // The table below was built from the components that were just replaced.
        if loaded {
            ui.ctx().request_repaint();
            return;
        }
        let lib = &lib;
        ui.horizontal(|ui| {
            ui.label("Find substance");
            ui.add(
                egui::TextEdit::singleline(&mut self.filter)
                    .hint_text("Name in the table and the library"),
            );
            if !self.filter.is_empty() && ui.button("Clear").clicked() {
                self.filter.clear();
                self.looked_up = None;
            }
        });
        let want = self.filter.to_lowercase();
        if !want.is_empty() {
            self.lookup(ui, lib, &want);
        }
        let selected = self.found_list(ui, run, lib);
        if !self.note.is_empty() {
            ui.label(&self.note);
        }

        let mut order: Vec<usize> = (0..rows.len())
            .filter(|i| want.is_empty() || rows[*i].name.to_lowercase().contains(&want))
            .collect();
        if let Some((col, desc)) = self.sort {
            order.sort_by(|a, b| {
                let o = compare(&rows[*a], &rows[*b], col);
                if desc { o.reverse() } else { o }
            });
        }
        ui.label(format!("{} of {} components", order.len(), rows.len()));
        if order.is_empty() && !rows.is_empty() {
            ui.weak("No component matches. Clear the search to see all.");
        }
        let clicked = table(ui, &mut self.sort, &rows, &order, run, lib).or(selected);
        if let Some(k) = clicked {
            run.component = Some(k);
            run.scan = run.components[k].apex_scan as usize;
        }
        self.cache = Some((stamp, rows));
    }

    /// Two buttons: save the components of the run as a `.cmp` file the old program reads, or load such
    /// a file in place of them. Returns whether components were loaded.
    fn saved_buttons(&mut self, ui: &mut egui::Ui, run: &mut Run) -> bool {
        if ui
            .button("Save components…")
            .on_hover_text("Writes the components and their peaks as a .cmp file (old format)")
            .clicked()
            && let Some(p) = dialog(run).save_file()
        {
            self.note = match save_components(run, &p) {
                Ok(n) => format!("Saved {n} components to {}", p.display()),
                Err(e) => format!("Save failed: {e}"),
            };
        }
        let busy = run.analysis.busy();
        if ui
            .add_enabled(!busy, egui::Button::new("Load components…"))
            .on_hover_text("Replaces the components with those of a .cmp file of this run")
            .on_disabled_hover_text("Wait until the analysis has finished")
            .clicked()
            && let Some(p) = dialog(run).pick_file()
        {
            let result = load_components(run, &p);
            self.note = match &result {
                Ok(n) => format!("Loaded {n} components from {}", p.display()),
                Err(e) => format!("Load failed: {e}"),
            };
            return result.is_ok();
        }
        false
    }

    /// Library names that contain the search text; a click lists the components that have it as a hit.
    fn lookup(&mut self, ui: &mut egui::Ui, lib: &Library, want: &str) {
        let names = lib
            .compounds
            .iter()
            .enumerate()
            .filter(|(_, c)| c.name.to_lowercase().contains(want))
            .take(8);
        ui.horizontal_wrapped(|ui| {
            ui.label("Library:");
            for (i, c) in names {
                if ui
                    .selectable_label(self.looked_up == Some(i), &c.name)
                    .clicked()
                {
                    self.looked_up = Some(i);
                }
            }
        });
    }

    /// Components with a hit for the looked-up substance, best first; returns the clicked one.
    fn found_list(&self, ui: &mut egui::Ui, run: &Run, lib: &Library) -> Option<usize> {
        let (i, ((..), id, _)) = (self.looked_up?, run.ident.as_ref()?);
        let found = find_substance(id, &lib.compounds, i);
        let mut clicked = None;
        ui.horizontal_wrapped(|ui| {
            ui.strong(&lib.compounds[i].name);
            if found.is_empty() {
                ui.label("not found in this run");
            }
            for (k, score) in found.iter().take(12) {
                let at = run.data.scan_times_s[run.components[*k].apex_scan as usize] / 60.0;
                if ui.link(format!("{at:.2} min: {score:.3}")).clicked() {
                    clicked = Some(*k);
                }
            }
        });
        clicked
    }
}

/// Rows of the active run; the RI comes from the calibration the identification was made with.
fn build_rows(run: &Run) -> Vec<ResultRow> {
    let Some(((lib, ..), id, shifted)) = &run.ident else {
        return Vec::new();
    };
    let times = &run.data.scan_times_s;
    let ris = (shifted.as_ref())
        .or(run.ident.as_ref().and_then(|i| i.0.1.as_deref()))
        .map(|c| c.component_ris(times, &run.components));
    result_rows(
        &run.components,
        times,
        ris.as_deref(),
        &run.peak_masses,
        id,
        &lib.compounds,
    )
}

/// Click order of a column header: ascending, descending, off.
fn next_sort(sort: Option<(usize, bool)>, col: usize) -> Option<(usize, bool)> {
    match sort {
        Some((s, false)) if s == col => Some((col, true)),
        Some((s, true)) if s == col => None,
        _ => Some((col, false)),
    }
}

/// The table of `order` (indices into `rows`); returns the component clicked in it.
fn table(
    ui: &mut egui::Ui,
    sort: &mut Option<(usize, bool)>,
    rows: &[ResultRow],
    order: &[usize],
    run: &Run,
    lib: &Library,
) -> Option<usize> {
    let mut clicked = None;
    ui.scope(|ui| {
        ui.style_mut().wrap_mode = Some(egui::TextWrapMode::Extend);
        TableBuilder::new(ui)
            .id_salt("results table")
            .striped(true)
            .sense(egui::Sense::click())
            .column(Column::auto().at_least(180.0))
            .columns(Column::auto(), 6)
            .column(Column::remainder())
            .header(20.0, |mut h| {
                for (c, title) in HEAD.iter().enumerate() {
                    h.col(|ui| {
                        ui.with_layout(align(c), |ui| header_cell(ui, sort, c, title));
                    });
                }
            })
            .body(|body| {
                body.rows(18.0, order.len(), |mut row| {
                    let r = &rows[order[row.index()]];
                    row.set_selected(run.component == Some(r.component));
                    for (c, text) in cells(r).into_iter().enumerate() {
                        row.col(|ui| {
                            ui.with_layout(align(c), |ui| {
                                // An unmatched component reads as a placeholder.
                                match r.compound.is_none() && c == 0 {
                                    true => ui.weak(text),
                                    false => ui.label(text),
                                };
                            });
                        });
                    }
                    let resp = row.response();
                    if resp.clicked() {
                        clicked = Some(r.component);
                    }
                    resp.on_hover_ui(|ui| hover(ui, r, run, lib));
                });
            });
    });
    clicked
}

fn header_cell(ui: &mut egui::Ui, sort: &mut Option<(usize, bool)>, col: usize, title: &str) {
    let arrow = match *sort {
        Some((s, d)) if s == col => [" ▲", " ▼"][d as usize],
        _ => "",
    };
    let t = egui::RichText::new(format!("{title}{arrow}")).strong();
    let resp = ui
        .add(egui::Label::new(t).sense(egui::Sense::click()))
        .on_hover_cursor(egui::CursorIcon::PointingHand)
        .on_hover_text("Click to sort: ascending, descending, off");
    if resp.clicked() {
        *sort = next_sort(*sort, col);
    }
}

/// Numbers (RT to RI score) are right-aligned so digits line up; name and ions stay left.
fn align(col: usize) -> egui::Layout {
    match (1..=6).contains(&col) {
        true => egui::Layout::right_to_left(egui::Align::Center),
        false => egui::Layout::left_to_right(egui::Align::Center),
    }
}

/// Cell texts of one row, as the old table showed them (a missing score as a dash).
fn cells(r: &ResultRow) -> [String; 8] {
    let num = |v: Option<f64>, digits: usize| v.map_or("–".into(), |v| format!("{v:.digits$}"));
    let ions: Vec<String> = r.quant_ions.iter().map(f32::to_string).collect();
    [
        r.name.clone(),
        format!("{:.2}", r.rt_s / 60.0),
        num(r.ri, 1),
        num(r.area, 0),
        num(r.overall, 3),
        num(r.spec_score, 3),
        num(r.ri_score, 3),
        ions.join(" "),
    ]
}

/// Missing values sort first.
fn compare(a: &ResultRow, b: &ResultRow, col: usize) -> std::cmp::Ordering {
    let num = |v: Option<f64>| v.unwrap_or(f64::NEG_INFINITY);
    match col {
        0 => a.name.to_lowercase().cmp(&b.name.to_lowercase()),
        1 => a.rt_s.total_cmp(&b.rt_s),
        2 => num(a.ri).total_cmp(&num(b.ri)),
        3 => num(a.area).total_cmp(&num(b.area)),
        4 => num(a.overall).total_cmp(&num(b.overall)),
        5 => num(a.spec_score).total_cmp(&num(b.spec_score)),
        6 => num(a.ri_score).total_cmp(&num(b.ri_score)),
        _ => a.quant_ions.len().cmp(&b.quant_ions.len()),
    }
}

/// Old tooltip: scores, RI, RT and the alternative hits of the component.
fn hover(ui: &mut egui::Ui, r: &ResultRow, run: &Run, lib: &Library) {
    let Some(((..), id, _)) = &run.ident else {
        return;
    };
    ui.strong(&r.name);
    ui.label(format!("RT {:.2} min", r.rt_s / 60.0));
    if let Some(ri) = r.ri {
        ui.label(format!("RI {ri:.1}"));
    }
    for h in id.hits[r.component].iter().take(5) {
        let mark = if Some(h.compound) == r.compound {
            "→ "
        } else {
            "   "
        };
        ui.label(format!(
            "{mark}{:.3}  {}",
            h.overall, lib.compounds[h.compound].name
        ));
    }
}

/// The file dialog for the saved components of `run`: `.cmp` next to the run, where it is found again.
fn dialog(run: &Run) -> rfd::FileDialog {
    let d = rfd::FileDialog::new()
        .add_filter("Saved components", &["cmp"])
        .set_file_name(format!("{}.cmp", run.name));
    match run.path.parent() {
        Some(dir) => d.set_directory(dir),
        None => d,
    }
}

/// Writes the components of `run` to `path`; returns how many.
fn save_components(run: &Run, path: &std::path::Path) -> Result<usize, String> {
    // The RI is that of the identification, so the old program lists it too.
    let mut ris = vec![-1.0; run.components.len()];
    for r in build_rows(run) {
        ris[r.component] = r.ri.unwrap_or(-1.0);
    }
    let compounds = md_core::cmp::compounds_of_components(
        &run.components,
        &run.peak_masses,
        &run.peaks,
        &run.data.scan_times_s,
        &ris,
    );
    md_core::cmp::write_cmp(path, &compounds).map_err(|e| e.to_string())?;
    Ok(compounds.len())
}

/// Replaces the components of `run` by those in `path`, which must be of this run; returns how many.
/// The selection and hits belong to the old components and go.
fn load_components(run: &mut Run, path: &std::path::Path) -> Result<usize, String> {
    let (components, masses, peaks) =
        md_core::cmp::components_from_file(path, &run.data.scan_times_s)
            .map_err(|e| e.to_string())?;
    let n = components.len();
    (run.components, run.peak_masses, run.peaks) = (components, masses, peaks);
    (run.component, run.ident, run.ident_job) = (None, None, None);
    run.rev += 1;
    Ok(n)
}

/// Asks for a file and writes the table; returns the note to show.
fn export(run: &str, rows: &[ResultRow]) -> String {
    let Some(p) = rfd::FileDialog::new()
        .add_filter("CSV", &["csv"])
        .set_file_name(format!("{run}-results.csv"))
        .save_file()
    else {
        return String::new();
    };
    match std::fs::File::create(&p)
        .map(std::io::BufWriter::new)
        .and_then(|f| write_csv(f, rows))
    {
        Ok(()) => format!("Saved {}", p.display()),
        Err(e) => format!("Export failed: {e}"),
    }
}

#[cfg(test)]
pub(crate) mod tests {
    use std::sync::Arc;

    use md_core::Chromatogram;
    use md_core::deconv::{Component, FittedPeak};
    use md_core::library::LibraryCompound;
    use md_core::search::{Identification, LibraryHit, PreparedLibrary, SearchSettings};

    use super::*;
    use crate::settings;

    pub(crate) fn run() -> Run {
        let s = SearchSettings::default();
        let compounds = vec![LibraryCompound::new("Alanine")];
        let lib = Arc::new(Library {
            path: "lib.mdlib".into(),
            name: "lib.mdlib".into(),
            prepared: PreparedLibrary::new(&compounds, &s),
            compounds,
            settings: s,
        });
        let hit = LibraryHit {
            compound: 0,
            spec_sim: 0.9,
            ri_sim: None,
            overall: 0.9,
        };
        let comp = |apex| Component {
            first_scan: 0,
            apex_scan: apex,
            last_scan: 2,
            model: Vec::new(),
            spectrum: Vec::new(),
            peaks: vec![FittedPeak {
                peak: 0,
                model_share: 1.0,
                integral: 5.0,
            }],
        };
        Run {
            path: "a.cdf".into(),
            name: "a".into(),
            data: Arc::new(Chromatogram {
                scan_times_s: vec![0.0, 60.0, 120.0],
                offsets: vec![0; 4],
                mz: Vec::new(),
                intensity: Vec::new(),
            }),
            minutes: vec![0.0, 1.0, 2.0],
            tic: vec![0.0; 3],
            sic: Vec::new(),
            sic_for: None,
            sic_job: None,
            components: vec![comp(1), comp(2)],
            peak_masses: vec![73],
            peaks: Vec::new(),
            saved_calibration: None,
            analysis: settings::Analysis::default(),
            scan: 0,
            component: None,
            visible: true,
            ident: Some((
                (lib, None, None),
                Identification {
                    hits: vec![vec![hit], Vec::new()],
                    assigned: vec![Some(hit), None],
                },
                None,
            )),
            ident_job: None,
            rev: 0,
        }
    }

    #[test]
    fn window_draws_every_row_and_caches_them() {
        let ctx = egui::Context::default();
        let mut runs = vec![run()];
        let mut r = Results {
            open: true,
            filter: "ala".into(),
            looked_up: Some(0),
            ..Default::default()
        };
        for _ in 0..3 {
            let mut out = ctx.run_ui(Default::default(), |ui| {
                r.window(ui.ctx(), &mut runs, 0);
            });
            out.textures_delta.clear();
        }
        let rows = &r.cache.as_ref().unwrap().1;
        assert_eq!(rows.len(), 2);
        assert_eq!(rows[0].name, "Alanine");
        assert_eq!(rows[0].area, Some(5.0));
    }

    #[test]
    fn saved_components_load_again_with_their_peaks() {
        use md_core::peaks::Peak;

        let mut a = run();
        for k in &mut a.components {
            k.model = vec![1.0, 4.0, 1.0];
        }
        a.peaks = vec![Peak {
            row: 3,
            first_scan: 0,
            apex_scan: 1,
            last_scan: 2,
            precise_apex: Some(1.25),
            signal_to_noise: 12.0,
            quality: 0.5,
            sharpness: 3.0,
        }];
        let path = std::env::temp_dir().join(format!("md-app-save-{}.cmp", std::process::id()));
        assert_eq!(save_components(&a, &path), Ok(2));

        let mut b = run();
        b.component = Some(1);
        let rev = b.rev;
        assert_eq!(load_components(&mut b, &path), Ok(2));
        std::fs::remove_file(&path).ok();
        // Both components shared one peak; the file gives each its own.
        let mut want = a.components.clone();
        want[1].peaks[0].peak = 1;
        assert_eq!(b.components, want);
        assert_eq!((b.peak_masses.clone(), b.peaks.len()), (vec![73, 73], 2));
        assert_eq!(
            (b.peaks[0].signal_to_noise, b.peaks[0].precise_apex),
            (12.0, Some(1.25))
        );
        assert!(b.component.is_none() && b.ident.is_none() && b.rev == rev + 1);

        // A file of another run is refused and leaves the run as it was.
        let mut other = run();
        other.data = std::sync::Arc::new(Chromatogram {
            scan_times_s: vec![0.0, 61.0, 120.0],
            ..(*other.data).clone()
        });
        let path = std::env::temp_dir().join(format!("md-app-save2-{}.cmp", std::process::id()));
        save_components(&a, &path).unwrap();
        assert!(load_components(&mut other, &path).is_err());
        std::fs::remove_file(&path).ok();
        assert_eq!(other.components.len(), 2);
    }

    #[test]
    fn a_new_result_rebuilds_the_rows() {
        let ctx = egui::Context::default();
        let mut runs = vec![run()];
        let mut r = Results {
            open: true,
            ..Default::default()
        };
        let draw = |r: &mut Results, runs: &mut Vec<Run>| {
            ctx.run_ui(Default::default(), |ui| r.window(ui.ctx(), runs, 0))
                .textures_delta
                .clear();
        };
        draw(&mut r, &mut runs);
        assert_eq!(r.cache.as_ref().unwrap().1.len(), 2);
        // The hits shrink with the components; a stale cache would index past them.
        runs[0].components.truncate(1);
        let id = &mut runs[0].ident.as_mut().unwrap().1;
        id.hits.truncate(1);
        id.assigned.truncate(1);
        runs[0].rev += 1;
        draw(&mut r, &mut runs);
        assert_eq!(r.cache.as_ref().unwrap().1.len(), 1);
    }
}
