# Contributors

Start with this file and .vendor/workflow-kit/AGENT_RULES.md.
Read the applicable sections of .vendor/workflow-kit/docs/CONTRIBUTING.md once when needed;
reuse them and reread only when relevant rules change. Do not recursively load
every linked document, generated skill or provider directory.
Reuse unchanged context/evidence; refresh live issue/dependency/ownership state.
The GitHub Project recorded in
`.github/workflow-project.json` owns order/status; issues own scope,
acceptance, dependencies and evidence; PRs own review. No second task ledger,
claim PRs or central steering issue.
Assign every issue a milestone, an explicit Project Priority and matching labels,
including Backlog and follow-ups; retain metadata after closure.
See .vendor/workflow-kit/docs/CONTRIBUTING.md#milestones for scope and the ongoing hotfix milestone.

New issues start in Backlog. Human-authorized triage may promote feasible issues
to Ready even while dependencies or external blockers remain. Ready is not proof
that blockers are resolved; resolve them before any start.

Start implementation only from Ready under the start policy in
.vendor/workflow-kit/docs/CONTRIBUTING.md#board-and-ownership (default: an explicit human request
for that work). If Ready is empty, analyze Backlog and propose the next executable
issue with rationale/blockers; do not implement or change status. A human must
authorize promotion; record readiness before moving to In progress (Doing).
Existing authorization covers review fixes within the same scope, not new tasks.

Before work, check issue, dependencies, current main and competing PRs. Use
codex/<issue>-<topic>, one driver and linked branch/PR. Assignment is not a lock.
Create issue-linked branches natively; verify PR/issue links after creation and
before review under .vendor/workflow-kit/docs/CONTRIBUTING.md#issue-branch-and-pr-links.
Branch names or mentions alone do not prove linkage; partial work must not close issues.
Preserve foreign work/ownership; no force-push or admin bypass. Human review/merge
remain default. Driver owns coordination, integration, proof and delivery;
delegation grants no ownership, scope expansion or approval authority.

Finish authorized scope through verification and PR delivery, not the entire
board. Before reporting agent work complete, await automatic reviews triggered by
Ready for Review in Automated review; fix or link actionable follow-ups for every
finding, then move to Human review for human acceptance. Return
PR to Draft/In progress for further work, then Ready for Review/Automated review
and await reviews again.
Confirmed review unavailability (e.g. exhausted tokens/quota) permits agent
completion only with the missing review and cause disclosed; pending/unknown is
not unavailable. See .vendor/workflow-kit/docs/CONTRIBUTING.md#publication-and-review. Prepare concrete
results before requesting missing approval; cite the blocking rule and continue
unaffected work.

Read architecture/contracts before behavior changes. Product/technology belongs
in design docs and milestone/issues. Accepted API changes update contracts and
consumers together; today's API is not immutable. Preserve validation,
authorization, privacy, integrity and attribution; keep secrets out of clients/logs.

Follow the shared engineering rules for the smallest complete change. Run focused
behavior and required issue/CI checks; never weaken acceptance or equate mocks with
live proof. Commit the kit's required managed outputs: provider skills, companions,
hook definitions and discovery, verified against pinned sources. Keep personal
snapshots, private artifacts, caches and replaced copies out of commits.
Record blocker/recovery and preserve
branch/PR. Human review requires completed automatic reviews (or a disclosed
confirmed outage), dispositioned findings and passed checks; Done requires
acceptance and merge. Disposition concrete findings with a fix or actionable linked
follow-up. Update durable docs; report results and uncertainty concisely.

## Hook preflight

At session start, check Ponytail and Impeccable hook status in the current harness
following [hook rules](.vendor/workflow-kit/AGENT_RULES.md#hooks). Warn once with activation steps
for confirmed failures or missing/disabled/untrusted hooks. If status is inaccessible,
continue without a warning or user check request; do not claim automatic execution.
Files or manual script runs do not prove live activation. Continue unaffected work.

Skill generation and hook installation belong to reviewed checkout provisioning
or explicit kit/tool upgrades, not routine task/turn preflight. Reuse a prepared
checkout; rerun affected validation after inputs change. Native hook trust stays
an explicit personal step; missing prerequisites do not authorize bypasses.

## Token and Actions budgets

Treat tokens and GitHub Actions usage as limited budgets, not free capacity.

- Read/search only relevant context; reuse valid evidence, keep output concise,
  and avoid duplicate tool calls, checks or unchanged-state polling.
- Prefer smaller, cheaper models when sufficient for a bounded task. Delegate
  only when expected token savings or useful parallel progress outweigh context,
  coordination and verification overhead; handle trivial work directly. Pass only
  necessary context, avoid duplicate work, and preserve requested model constraints.
- Prefer existing workflows or native GitHub features over new automation.
  Before adding workflows/jobs or expanding triggers, schedules, matrices, runners
  or artifact retention, record necessity and estimated incremental monthly usage
  and cost in the issue/PR, including run frequency and billing assumptions.
  Verify fit within the existing budget; unknown or insufficient budget requires
  explicit user approval before enabling additional usage. Convenience alone
  does not justify recurring cost; never raise spending limits autonomously.
- Scope runs safely, cancel superseded runs when safe, and reuse valid results.
  Preserve required checks, security and acceptance; distinguish estimates from
  measured savings. Do not add a monitoring workflow just to enforce this policy.

## Code Review Rules

Preserve validation, authorization, privacy, integrity and attribution. Fix root
causes at the shared boundary. Keep producer, consumer, tests and documentation
consistent. Reuse standard/native facilities before introducing dependencies.

## Project-specific contracts

Read local CONTRIBUTING.md alongside the shared rules. Record architecture,
build/test commands and acceptance in project docs before behavior changes.
