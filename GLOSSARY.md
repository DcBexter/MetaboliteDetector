# Glossar

Domänenbegriffe des Projekts, wie sie in Code, Issues und docs/architecture.md benutzt werden.
Architekturbegriffe (module, interface, seam, adapter, deep/shallow, leverage, locality)
folgen dem Skill `codebase-design` und werden nicht übersetzt.

| Begriff | Bedeutung | Code |
|---------|-----------|------|
| **Run** | Ein GC-MS-Lauf: eine Datei (ANDI-MS `.cdf`, FastFlight, alter Cache) mit Scans, Zeiten und Spektren. | `md_core::Chromatogram` |
| **Derivative cache** | Erste Ableitung von TIC/SIC, je Run einmal beim ersten Zeichnen berechnet; das SIC-Feld wird mit einem neuen SIC geleert; beide werden freigegeben, sobald die Ableitung aus ist. | `Run::d_tic`, `Run::d_sic` (md-app) |
| **Scan** | Ein Spektrum zu einem Zeitpunkt des Runs; m/z innerhalb des Scans sortiert. | `Chromatogram` (flache Arrays mit Scan-Offsets) |
| **Ion trace** | Intensität einer (nominalen oder exakten) Masse über alle Scans. Einmal pro Run als ion-major Index gebaut, dicht oder (bei dünn besetzten Exact-Mass-Runs) sparse. | `md_core::peaks::IonTraces`, `trace(row, scans, buf)` |
| **Peak** | Ein Maximum in einem Ion trace mit Grenzen, Apex, interpoliertem Apex und S/N. | `md_core::peaks::Peak` |
| **Peak pass** | Welche Peak-Erkennung ein Aufrufer will: *Reordered* (Hauptfenster und Batch, `detect_all_peaks` mit `reorder_peaks` auf Exact-Mass-Runs) oder *Fresh* (Rescan, `detect_peaks` wie das alte `ResultSeries::scanIons`). | `md_core::analysis::PeakPass` |
| **Component** | Ein deconvolvierter Stoff: Peaks mehrerer Ion traces mit gemeinsamer Retentionszeit, Modellpeak und Spektrum. | `md_core::deconv::Component` |
| **Analyze run** | Die Kette Ion traces, Peaks, Deconvolution für einen Run hinter einer Interface. | `md_core::analysis::analyze_run`; Peak-Stufe allein: `find_run_peaks` |
| **Library** | Sammlung von Referenzspektren mit RI, Quant-Ionen und FRAG-Angaben (`.lbr`, `.msl`, `.mdlib`). | `md_core::library` |
| **Library search** | Zuordnung von Components zu Library-Einträgen über Spektrum- und RI-Score, greedy eins zu eins (Heap über die besten Treffer je Component, keine globale Paarliste). | `md_core::search` |
| **Retention index (RI)** | Zeitnormierung über Alkan-Marker; Kalibrierung RT zu RI. | `md_core::ri` |
| **Study** | Viele Runs einer Fragestellung, die gemeinsam quantifiziert werden: Runs, Library, Scan- und Quant-Settings, Rescan ja/nein. Ein Aufruf für Batch-Fenster, MID-Fenster und md-cli. | `md_core::batch::{Study, run_study, StudyResult}` |
| **Batch** | Die Quantifizierung einer Study: Substanz × Run Tabelle mit Quant-Ionen. | `md_core::batch::{run_study, table}` |
| **Rescan** | Zweiter Pass des Batch: fehlende Substanzen werden je Run an der erwarteten RT mit frischen Peaks nachgesucht. Nutzt die Fresh-Peak-Liste aus Pass 1, nicht den Ion-Index; liest einen Run nur, wenn dort ein Kandidat liegt. | `md_core::batch::rescan` |
| **Quant ion** | Die Masse, deren Peakfläche eine Substanz im Batch quantifiziert. | `QuantSettings` |
| **Normalisierung** | Umrechnung der Batch-Tabelle (interner Standard, Gesamtsignal, Pool, z-Score). | `md_core::normalize` |
| **MID** | Mass isotopomer distribution: Verteilung der Markierung je Fragment, korrigiert um natürliche Isotope. | `md_core::mid` |
| **Job** | Hintergrundarbeit der App mit Progress und Cancel; die UI pollt den Zustand. | `md_app::jobs::Job`; `md_core::job::{Cancel, Progress}` |
| **View model** | Fertig sortierte, gefilterte und formatierte Zeilen einer Tabelle, aus der Core-Tabelle plus Ansichtsoptionen gebaut und per Generation gecacht. | `md_core::view` (`BatchViewCache`) |
| **Legacy numeric** | Benannte Rundungen und Rechenwege, die 3.3 bit für bit nachbilden (C-Trunkierung `c_int`, `sg_trunc`, `round_sig8`, glibc `pow`/`hypot`). | `md-core/src/legacy_numeric.rs` (crate-intern) |
| **Golden test** | Vergleich der neuen Kette gegen Ergebnisse des alten Binary 3.3 auf den Beispiel-Runs. | `crates/md-core/tests/golden_old_md.rs` |
| **Stage fixture** | Zwischenausgabe einer Stufe von 3.3 (Peaks, Regionen, Components mit Modellpeak) je Beispiel-Run als TSV, für den Vergleich Stufe für Stufe. | `crates/md-core/tests/fixtures/legacy/`, Test `legacy_stages.rs` |
| **Settings** | Die `AN_*`/`LS_*`-Parameter der Analyse; eine Quelle in md-core. | `md_core::settings::GcmsSettings`, `KEYS` (Key-Tabelle: Schlüssel, Anzeigetext, Bereich) |
