---
name: ponytail-help
description: >
  Quick reference for ponytail levels, skills and commands. One-shot display.
  Use for /ponytail-help, "ponytail help", "how do I use ponytail".
---

# Ponytail Help

Display this reference card when invoked. One-shot, do NOT change mode,
write flag files, or persist anything.

## Levels

| Level | Trigger | What change |
|-------|---------|-------------|
| **Lite** | `/ponytail lite` | Build what was asked, name the smaller option in one line. |
| **Full** | `/ponytail` | The smallest complete change, a check where the logic needs one, and a reply that names what was skipped and any risk. Default. |
| **Ultra** | `/ponytail ultra` | Also questions the request and pushes back before building. |

Level sticks until changed or session end.

## Skills

| Skill | Trigger | What it does |
|-------|---------|--------------|
| **ponytail** | `/ponytail` | Lazy mode itself: least new code, clear replies that name skipped work and risks. |
| **ponytail-review** | `/ponytail-review` | Quality review of a diff: bugs, security, load, missing tests, speed, what to cut. Each finding says what goes wrong and how to fix it. |
| **ponytail-audit** | `/ponytail-audit` | The same quality review for the whole repo, ranked. |
| **ponytail-debt** | `/ponytail-debt` | Harvest `shortcut:` comments into a tracked ledger. |
| **ponytail-gain** | `/ponytail-gain` | Measured-impact scoreboard: less code, less cost, more speed. |
| **ponytail-help** | `/ponytail-help` | This card. |

Select the repository skill in your agent's skill picker. Codex CLI supports
`$ponytail`; desktop picker syntax can vary. Claude Code uses the slash-command
forms above. Other hosts can load the named skill or its `SKILL.md` explicitly.

## Deactivate

Say "stop ponytail" or "normal mode". With hooks, `/ponytail off` also works;
resume with `/ponytail` (off enables the default level; an active mode reports its level).

## Checkout scope

This checkout bundles skills plus Node.js lifecycle hooks for Codex, Claude Code,
Copilot and Cursor. From a reviewed checkout, initialize them with
`node .vendor/workflow-kit/scripts/install-ponytail-hooks.mjs`, then enable the host's project trust/hooks
to activate the default mode at session start. Default priority: `PONYTAIL_DEFAULT_MODE`, then personal
Ponytail config `defaultMode`, then `full`. Supported defaults: off/lite/full/ultra.
`/ponytail default lite` saves a personal default; ordinary mode switches do not.
No global statusline or automatic plugin update is installed. Other hosts use
the skills directly. Repository workflow, validation and security rules remain authoritative.

## Update

Updates arrive through repository changes. See
[checkout usage and updates](https://github.com/snagg-app/workflow-kit/blob/main/README.md#install-or-update-in-a-project).
Regenerate local links and cloud files together; retain patches and MIT attribution.

## More

Full docs + examples: https://github.com/DietrichGebert/ponytail
