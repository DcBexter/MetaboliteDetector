# Sicherheit und automatische Wartung

Schwachstellen privat über **Security → Report a vulnerability** melden.
Keine Laborrohdaten, persönlichen Informationen oder Zugangsdaten veröffentlichen.

GitHub prüft Rust, JavaScript/TypeScript und Actions mit CodeQL; Dependabot meldet
bekannte Schwachstellen in erkannten Abhängigkeiten und erstellt Sicherheitsupdates.
Dependency Review ergänzt neue Abhängigkeiten im PR. Rust-Tests und Sonar bleiben
bestehende Pflichtprüfungen. Secret Scanning und Push Protection sind eingeschaltet.

Renovate verwaltet Cargo, Actions, den Rust-Toolchain-Pin und vorhandene Submodule.
Alle Update-Typen einschließlich Major-, Sicherheits- und Lockfile-Updates dürfen
nach grünen Pflichtchecks automatisch mergen. `ignoreTests` bleibt aus; der Branch
muss aktuell sein. Neue Versionen warten wie bisher drei Tage. Native GitHub-Checks
sind an ihre Apps gebunden. Diese Regel folgt dem ausdrücklichen Nutzerauftrag vom
24.10.2026; fachliche Validierung ist kein automatisch erfülltes Gate.

Die Tests ohne bereitgestellte Originaldaten ersetzen insbesondere nicht die
menschliche Labordaten-Abnahme in [#147](https://github.com/DcBexter/MetaboliteDetector/issues/147).
Keine Laborrohdaten werden für diese Automation beschafft, hochgeladen oder erzeugt.
Keine bezahlten KI-Funktionen oder neuen Release-/Deployment-Trigger.
