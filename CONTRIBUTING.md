# Contributing

Follow [the shared workflow](.vendor/workflow-kit/docs/CONTRIBUTING.md) and
[agent rules](.vendor/workflow-kit/AGENT_RULES.md). Put project-specific
architecture and build/test commands in `AGENTS.md` and the project documentation.

## This project

Layout, build and test commands: [docs/architecture.md](docs/architecture.md#layout).
New crates go under `crates/` and inherit versions, dependencies and lints from the workspace.

Code analysis runs in CI ([docs/architecture.md](docs/architecture.md#code-analysis)). A maintainer sets it up
once: switch off Automatic Analysis in SonarQube Cloud (Administration → Analysis Method) and add a
SonarQube Cloud token as the repository secret `SONAR_TOKEN`.
