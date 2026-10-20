# Session brief

What every agent session on this project needs to know. A session prompt names only the issue and what is
specific to it; the rest is here. Shared rules live in the workflow kit and are linked, not copied.

## Goal and where things are

- Goal: a Rust + egui rewrite of the old GC-MS tool with full parity, then better ([PRODUCT.md](../../PRODUCT.md),
  [feature inventory](../reverse-engineering/feature-inventory.md)). Performance on huge datasets is a priority.
- Layout, build and test commands: [architecture.md](../architecture.md#layout) and
  [Build and test](../architecture.md#build-and-test). `md-core` has no UI dependencies, `md-app` draws.
- Design rules for the UI: [DESIGN.md](../../DESIGN.md) and the [UI checklist](../ui-quality.md).

## Rules for this project

- No example data and nothing from `Old/` or generated old-version files go into a commit (both are in `.gitignore`).
  New downloads need the user's permission.
- Issues, PRs and handoff comments are written in simple German ([workflow kit](../../.vendor/workflow-kit/docs/CONTRIBUTING.md#einfache-sprache)).
  Code, code comments and `docs/` are English.
- The PR body keeps `Closes #N` and the line `Agent: claude, Session: ID` that `board.mjs start` writes
  ([PR body](../../.vendor/workflow-kit/docs/CONTRIBUTING.md#delivery)). Commits and PR bodies end with the attribution lines of the harness.
- Do not change settings: no GitHub, SonarCloud or ruleset changes, no tokens, no force-push. A needed change becomes an issue.
- Merge `origin/main` into your branch often, and before the final push.

## Flow

1. `node .vendor/workflow-kit/scripts/board.mjs start ISSUE` (only for issues in Ready), then work on the branch it opens.
   In a fresh worktree run `scripts/link-example-data.ps1` (or `.sh`) before testing: without the example data the golden tests
   only print "skipped" and stay green.
2. `cargo fmt --all`, `cargo clippy --workspace --all-targets -- -D warnings` and `cargo test --workspace` are green. A UI change also
   gets the screenshot check from the [UI checklist](../ui-quality.md). It runs locally after `scripts/link-example-data.ps1` (or `.sh`):
   `cargo test --release -p metabolite-detector screenshots -- --ignored` (about 9 s, 42 images). Skip it only with the real error message in the handoff.
3. `board.mjs done ISSUE PR` waits for CI and the automatic review; then `board.mjs done ISSUE PR FILE` with the handoff comment, then `board.mjs merge PR`.
4. Open Sonar findings: `scripts/sonar-open.sh [PR]`.
5. If CodeRabbit is rate limited, start a fresh `studio-core:reviewer` agent on the fixed head and say so in the handoff
   and say so in the handoff. Fix every finding or link a follow-up issue: every unfixed finding gets its own follow-up
   (`board.mjs new`), minor ones too. Without a replacement reviewer, the handoff says "Eingeschränkte Reviewer" instead.
6. The GraphQL quota is shared by all sessions. When it runs out, use `board.mjs quota-wait` instead of retrying.

## Hotspots: one driver at a time

These files are changed by many sessions at once, and most merge conflicts happen here. Before editing one, look at the open PRs
(`gh pr list`) and the linked sessions. If another driver has it, wait for that PR to merge or ask that session.

| File | Why |
|---|---|
| `crates/md-app/src/central.rs` | the main window: toolbar, plots, panels |
| `crates/md-app/src/layout.rs`, `style.rs` | window placement, theme, fonts |
| `crates/md-app/src/settings.rs` | every analysis setting |
| `.github/workflows/ci.yml` | the CI jobs |

Cross-cutting UI work (theme, layout, spacing across windows) belongs to the UI-quality session. A session that owns one window only
changes that window and tells the UI-quality session about shared needs.
