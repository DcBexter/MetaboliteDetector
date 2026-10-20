# Contributing

Follow [the shared workflow](.vendor/workflow-kit/docs/CONTRIBUTING.md) and
[agent rules](.vendor/workflow-kit/AGENT_RULES.md). Put project-specific
architecture and build/test commands in `AGENTS.md` and the project documentation.

Agent sessions: [session brief and hotspot files](docs/agents/session-brief.md).

## This project

Layout, build and test commands: [docs/architecture.md](docs/architecture.md#layout).
New crates go under `crates/` and inherit versions, dependencies and lints from the workspace.

Code analysis runs in CI ([docs/architecture.md](docs/architecture.md#code-analysis)). A maintainer sets it up
once: switch off Automatic Analysis in SonarQube Cloud (Administration → Analysis Method) and add a
SonarQube Cloud token as the repository secret `SONAR_TOKEN`.

Open Sonar issues of main, or of pull request N: `scripts/sonar-open.sh [N]` (needs Node, no token).

Release: raise `version` in the root `Cargo.toml` on main, then tag that commit `vX.Y.Z` and push the tag.
`.github/workflows/release.yml` builds the portable Windows exe, the MSI (`installer/metabolite-detector.wxs`,
WiX 6) and the Linux tar.gz, and publishes them as a pre-release. MSI versions are numeric, so `vX.Y.Z-rc1`
installs as `X.Y.Z`.
