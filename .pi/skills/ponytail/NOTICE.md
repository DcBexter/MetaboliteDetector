# Ponytail provenance

The six sibling `ponytail*` skill packages and shared `.agents/hooks/ponytail-*.js` derive from
[DietrichGebert/ponytail 5.2.0](https://github.com/DietrichGebert/ponytail/tree/v5.2.0),
the pinned `.vendor/ponytail` submodule.
Copyright (c) 2026 DietrichGebert; [MIT license](LICENSE.md) applies to all skills and hooks.

Vaultdex adaptations: remove the core's host-specific `argument-hint` metadata;
reuse required repository tests; make help accurate for checkout discovery and
hooks; keep debt findings in the issue workflow instead of a ledger file. Other
content is unchanged apart from line-ending normalization. All six host skill
copies must match.

Hook adaptations: native project manifests run an explicitly installed user-local
snapshot through a launcher that starts a fixed Node, and pass the host explicitly.
Missing snapshots emit setup guidance without executing checkout code. Runtime
state is separated by Git checkout and host in the user's Ponytail config
directory instead of sharing another installation's flag. The activation hook
omits global statusline setup and personal Claude-settings reads. Other hook
adaptations preserve live mode on resume/compact, deliver rules on Codex mode
switches and Copilot subagent creation, and accept Copilot's `agentName`. Stdin
fallbacks stay referenced until EOF or timeout, then release stdin and drain
stdout. Explicit hosts ignore inherited plugin/project variables. Other hook
logic is upstream; no npm dependency or runtime download is added.
