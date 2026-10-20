---
name: spec-review
description: "Check a spec issue against the code: complete, correct, as planned. Open decisions and architecture candidates are settled with the user through the question tool; agreed work becomes new sub-issues of the spec."
disable-model-invocation: true
argument-hint: "<spec issue number, URL or title words>"
---

# Spec review

Target spec: `$ARGUMENTS`. A number or URL names it directly; title words go through
`gh issue list --label spec --state all --search "<words>"`. Empty or several matches: ask once which spec.
Only a human starts this skill. Keep a list of every sub-issue created or newly linked during this invocation (steps 5 and 8), even if it closes before the run ends.

Every acceptance point and every planned decision of the spec ends with one **verdict**:

- **done**: implemented as planned, proof `file:line` on the review ref.
- **deviates**: implemented, but differently from the plan.
- **open**: missing or only partly there.

Project rules win: the project's `AGENTS.md` names its tracker, architecture and spec acceptance docs; read the rows
that apply before step 1. `board.mjs` is the kit's `scripts/board.mjs`.

## Steps

1. **Read the plan.** Spec body (to a file when large), comments, sub-issues with state
   (`gh api --paginate "repos/{owner}/{repo}/issues/SPEC/sub_issues?per_page=100" --jq '.[]|[.number,.state,.title]'`) and their merged PRs.
   The review ref is the remote branch those PRs merged into, after `git fetch`; note its commit.
   Done when every acceptance point and planned decision is listed.

2. **Architecture pass.** Read [`improve-codebase-architecture`](../improve-codebase-architecture/SKILL.md) and run its
   steps 1 and 2 (explore, HTML report) with the spec's code as the named direction, on the review ref.
   It is user-invoked, so follow the file; skip its closing question, step 7 asks about every candidate.
   A candidate counts as **deviates** when it breaks the spec, an ADR or the project's target architecture.
   Done when every candidate is listed for step 7.

3. **Verdict per point.** Check each listed point against the review ref: presence with `file:line`,
   absence with one `git grep` per pattern over the whole result. Done when every point carries a verdict and its proof.

4. **Collect the questions.** One question per **deviates** or **open** verdict whose fix needs a decision
   (two reasonable ways, where a rule lives, scope, accept the deviation or fix it), and one per architecture
   candidate from step 2, always: also when it breaks the spec, an ADR or the target architecture, and whatever
   strength you gave it. A gap with one obvious fix needs no question: it becomes a ticket in step 5.

5. **Gap tickets.** Each **open** gap with one obvious fix, without a matching open sub-issue, first gets a search of
   open issues (`gh issue list --search "<finding words>"`): an unlinked match from an earlier run is linked with
   `board.mjs sub SPEC N`. Otherwise it becomes a new issue via `board.mjs new` (required fields per the project's
   tracker doc), then `board.mjs sub SPEC NEW`. Right before each `board.mjs new`, repeat the duplicate check, because
   parallel runs and other agents file issues in the meantime: list the issues updated in the last 24 hours, open and
   closed (`gh api "repos/{owner}/{repo}/issues?state=all&since=<UTC now minus 24 h, ISO 8601>&per_page=100" --jq '.[]|select(.pull_request|not)|[.number,.state,.title]'`),
   and match them by the title words of the new ticket and by a reference to the spec in the body (`#SPEC` or the spec URL; read the body only of title candidates).
   A match, open or closed, is linked with `board.mjs sub SPEC N` instead of creating the issue. When `new` refuses because an open issue has the same title, link that
   issue instead (`board.mjs sub SPEC N`). No ticket carries an open decision; `needs-human-input` is only for
   work a human must do. Architecture candidates never get a ticket here.
   A ticket or acceptance point about a provider interface names the capability in general terms ("a source may create
   variants", an adapter contract capability) and the provider only as the first source or adapter, never as the subject.

6. **Report on the spec.** The branch and the sub-issues move during a run (other agents merge and file). Before
   writing, `git fetch` again and re-read the sub-issues: when the review ref moved past the noted commit or a
   sub-issue appeared, re-check every point whose files `git diff --name-only OLD NEW` lists and every new sub-issue,
   then note the new commit. The report names the commit it checked.
   One comment on the spec in plain language, written to a file and posted with
   `gh issue comment SPEC --body-file FILE` (also the addendum in step 8): verdict table with proofs, gap tickets from
   step 5, path of the HTML report, and the questions from step 4 that step 7 settles.

7. **Walk the user through.** Open the HTML report in the user's external browser (`start` on Windows, `open` on
   macOS, `xdg-open` on Linux) and give its path, so the user sees each candidate while answering; then call the Skill tool with "grilling";
   the questions from step 4 are its first frontier, in plain language (kit `docs/CONTRIBUTING.md`, Einfache Sprache).
   Here every round and the final confirmation go through the host's question tool (Claude Code: `AskUserQuestion`)
   instead of grilling's ❓/➡️ text: at most 4 questions per call, a bigger round takes further calls, your
   recommended option first and marked as recommended. The ❓/➡️ text is at most a short overview before the call. Each architecture candidate gets
   the options "Sub-issue now", "Sub-issue after the release", "Discuss first (grilling)" and "Report only", worded in
   the user's language; "Discuss first" opens its own branch of the tree. Every answer ends in exactly one outcome:
   follow-up work (its scope and timing), accepted deviation, rejected candidate, or left on purpose. A rejection with
   a load-bearing reason gets the ADR offer of `improve-codebase-architecture` step 3. Done when the frontier is empty
   and the user confirmed.

8. **Agreed work to tickets.** Right after step 7, repeat the freshness check of step 6 and narrow or drop agreed
   work that the new commits or sub-issues already did. Then each answer that asks for follow-up work gets its ticket the way
   step 5 does, including its duplicate check right before each `board.mjs new` (last 24 hours, title words and spec
   reference; a match is linked with `board.mjs sub SPEC N`); the body states the decision with date and source, and "after the release" goes into its target
   release. Adjust the spec's text and acceptance point to the decisions (`board.mjs body`) and add a short addendum
   comment (Nachtrag) on the spec: decisions, outcomes, every sub-issue created or newly linked in this invocation.
   Done when every agreed piece of work has an open sub-issue, read back from the sub-issues API.
   Human review requires all sub-issues closed, no gap left, and the invocation's created-or-linked list empty.
   Then prove acceptance in the project's spec handoff comment and run `board.mjs field SPEC Status "Human review" --spec-review`;
   the command re-reads the children immediately before writing. Otherwise keep or return the spec to Backlog
   (`board.mjs field SPEC Status Backlog --spec-review`) and name the new or linked sub-issues in the Nachtrag.
   It stays open: only a human closes it or sets Done (`AGENT_RULES.md`, Hard rules).
