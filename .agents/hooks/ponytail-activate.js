#!/usr/bin/env node
// ponytail — Claude Code SessionStart activation hook (also Codex, Copilot,
// Grok, CodeBuddy and Cursor sessionStart)
//
// Runs on every session start:
//   1. Writes the active mode to the runtime's user-local state directory
//   2. Emits ponytail ruleset as hidden SessionStart context

const { getDefaultMode } = require('./ponytail-config');
const { getPonytailInstructions } = require('./ponytail-instructions');
const {
  clearMode,
  cursorRuleNotice,
  cursorRulePath,
  isCodeBuddy,
  isCodex,
  isCopilot,
  isCursor,
  readMode,
  setMode,
  writeHookOutput,
} = require('./ponytail-runtime');


// Native matcher groups pass this flag for resume/compact (and Claude forks).
// No stdin read: activation must also work when the host leaves its pipe open.
const mode = process.argv[3] === 'continue' ? readMode() || 'off' : getDefaultMode();

// "off" mode — skip activation entirely, don't write flag or emit rules
if (mode === 'off') {
  clearMode();
  return;
}

// Cursor with the always-on rule in the workspace: the rule already carries the
// ruleset and would contradict any other level, so leave the flag alone and
// hand the model a one-line notice instead of a second copy (#817).
if (isCursor) {
  const rule = cursorRulePath();
  if (rule) {
    try {
      writeHookOutput('SessionStart', mode, cursorRuleNotice(rule));
    } catch (e) {
      // Silent fail — stdout closed/EPIPE at hook exit must not surface as a hook failure
    }
    return;
  }
}

// 1. Write flag file
try {
  setMode(mode);
} catch (e) {
  // Silent fail -- flag is best-effort, don't block the hook
}

// 2. Emit the ponytail ruleset, filtered to the active intensity level.
const output = getPonytailInstructions(mode);

// Vaultdex: no global statusline setup, personal-settings reads or codebase map (the kit packages no ponytail-map.js).

try {
  writeHookOutput('SessionStart', mode, output);
} catch (e) {
  // Silent fail — stdout closed/EPIPE at hook exit must not surface as a hook failure
}
