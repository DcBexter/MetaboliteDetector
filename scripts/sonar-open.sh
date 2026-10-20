#!/usr/bin/env sh
# Open SonarQube Cloud issues of main, or of pull request N: `scripts/sonar-open.sh [N]`.
# Prints the count, then one line per issue: `rule file:line message`. The project is public, so no token.
# Needs Node 18+ (the workflow already uses it) and works in Git Bash and on Linux.
set -eu
exec node -e '
const key = "DcBexter_MetaboliteDetector";
const pr = process.argv[1];
if (pr !== undefined && !/^[0-9]+$/.test(pr)) { console.error("usage: sonar-open.sh [PR number]"); process.exit(2); }
(async () => {
  const issues = [];
  for (let page = 1; ; page++) {
    const q = new URLSearchParams({ componentKeys: key, resolved: "false", ps: "500", p: String(page) });
    if (pr) q.set("pullRequest", pr);
    const res = await fetch("https://sonarcloud.io/api/issues/search?" + q);
    if (!res.ok) { console.error("SonarCloud answered " + res.status + ": " + (await res.text())); process.exit(1); }
    const body = await res.json();
    issues.push(...body.issues);
    if (issues.length >= body.total || body.issues.length === 0) break;
  }
  console.log(issues.length);
  for (const i of issues) console.log(i.rule + " " + i.component.replace(key + ":", "") + ":" + (i.line ?? 0) + " " + i.message);
})().catch((e) => { console.error(e.message); process.exit(1); });
' -- "$@"
