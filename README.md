# MetaboliteDetector

Rust rewrite of MetaboliteDetector, a tool for deconvolution, identification and quantification
of GC/MS metabolomics data. Runs natively on Windows, Linux and macOS without R or Qt installs.

- Architecture, performance rules, build and test commands: [docs/architecture.md](docs/architecture.md)
- Milestones: [docs/roadmap.md](docs/roadmap.md)
- What the old 3.3 version did (reverse-engineered): [docs/reverse-engineering/](docs/reverse-engineering/)

Example data (15 public GC-MS runs, ~90 MB, not in git): run `example-data/fetch.sh`, or
`.\example-data\fetch.ps1` on Windows. Then:

```sh
cargo run --release -p metabolite-detector -- example-data/TargetSearchData/inst/gc-ms-data/7235eg04.cdf
```

Pass several files (or pick several in File → Open) to overlay them; TIC, SIC (m/z ± tolerance) or both.

---

## Project template (workflow kit)

Use **Use this template** on GitHub to start a project with the shared
[Workflow Kit](https://github.com/snagg-app/workflow-kit). The kit and its Ponytail /
Impeccable sources stay pinned Git submodules; improvements arrive as update PRs.
Add your own application and project-specific contracts without copying kit code.

For agent-guided setup, give your agent the [setup entry](https://github.com/snagg-app/workflow-kit/blob/main/SETUP.md).
It asks for repository, board and integration choices, then performs and verifies setup.

## Setup

Use Git, Node 26 and authenticated GitHub CLI:

```sh
git clone --recurse-submodules https://github.com/OWNER/REPO.git
cd REPO
node .vendor/workflow-kit/scripts/init-project.mjs --existing
node .vendor/workflow-kit/scripts/setup-skills.mjs
node .vendor/workflow-kit/scripts/setup-github.mjs OWNER/REPO
```

All shared harness commands run from this project's root and act on its checkout.
This template contains no forwarding
scripts; keep only actual product-specific scripts when adding an application.

The final command copies the public board template, links this repository and
creates missing labels. Commit `.github/workflow-project.json`. Configure native
statuses in order: Backlog → Ready → In progress → Automated review → Human review → Done.
Ready for Review starts Automated review; Human review requires completed automatic
reviews, dispositioned findings and passed checks (disclose confirmed unavailable
reviewers under the shared workflow exception). Further edits restart the review cycle.
Done still requires human acceptance and merge. Configure native
Project Auto-add for your repo and verify status automations in the GitHub UI;
cancelled work must not become Done. No private Vaultdex items are copied.

After reviewing the checkout, install hooks explicitly:

```sh
node .vendor/workflow-kit/scripts/install-ponytail-hooks.mjs
node .vendor/workflow-kit/scripts/install-impeccable-hooks.mjs
```

Review/enable hook definitions in each agent. Personal settings and plugin installs
are preserved. CodeRabbit requires GitHub App access; Codex requires a connected
repository and automatic reviews enabled in Codex settings. Organization-wide App
access can already cover the new repo. Configuration alone does not prove activation.
Read the pinned kit's setup/security notes first; the known Ponytail PATH finding
vaultdex/workflow-kit#3 is not fixed by moving scripts into the kit.

## Working and updating

Follow [CONTRIBUTING](CONTRIBUTING.md) and [AGENTS.md](AGENTS.md). Fill in your
architecture, build/test commands and product acceptance. Product CI and Sonar
project keys are specific to your project; this template does not invent them.

Dependabot proposes kit/Actions updates weekly. Renovate is explicitly disabled to
avoid duplicate proposals in organizations where its App covers all repos. If you
prefer Renovate, replace the Dependabot config with the kit's documented configuration.

For an update PR, inspect the new source revision, then run:

```sh
git -C .vendor/workflow-kit fetch origin
git -C .vendor/workflow-kit checkout <reviewed-kit-commit>
git add .vendor/workflow-kit
git submodule update --init --recursive
node .vendor/workflow-kit/scripts/init-project.mjs . --existing
node .vendor/workflow-kit/scripts/setup-skills.mjs .
```

Commit generated provider skills, companions and managed configuration in the same
PR. Setup moves replaced skill copies into ignored `.workflow-kit/replaced/`;
init preserves foreign hook handlers. Personal hook installation and trust remain
separate steps, including renewed review of changed hook definitions.
CI regenerates managed configuration and discovery, stages new paths with
`git add --intent-to-add --all`, then checks `git diff --exit-code HEAD` so missing
or stale committed output fails. Stage the new gitlink before the submodule
update; otherwise Git checks out the old pin again.

For upgrades from a wrapper-based template, update any external workspace bootstrap
and hook-install caller before removing its old entrypoint. `init-project` removes
only unchanged owned wrappers and migrates recorded hook metadata; it preserves
foreign handlers. For Snagg Symphony, deploy vaultdex/snagg-symphony#10 first.
The kit's separate review-policy PR #24 is not included in this pin; apply it only
after its required human acceptance and merge.

The included PR check verifies kit integration, not application correctness.
Budget: one standard Linux job, no schedule/cache/artifacts, maximum five minutes
per run. Public standard-runner usage is free; review private-repo and reviewer
budgets before enabling new projects. No paid plans are enabled by setup.
