#!/usr/bin/env node
// Links the main checkout's example data into this worktree (junction on Windows, symlink elsewhere),
// so golden tests do not skip. Same as link-example-data.ps1/.sh, but runs where only node is allowed.
import { execFileSync } from 'node:child_process';
import { existsSync, lstatSync, symlinkSync } from 'node:fs';
import { dirname, join } from 'node:path';

const git = (...args) => execFileSync('git', args, { encoding: 'utf8' }).trim();
const root = git('rev-parse', '--show-toplevel');
const main = dirname(git('rev-parse', '--path-format=absolute', '--git-common-dir'));
const fetchHint = process.platform === 'win32' ? 'example-data/fetch.ps1' : 'example-data/fetch.sh';

for (const d of ['TargetSearchData', 'old-results']) {
  const dst = join(root, 'example-data', d);
  const src = join(main, 'example-data', d);
  // lstat also sees a dangling link, which existsSync does not
  if (existsSync(dst) || safeLstat(dst)) { console.log(`ok: ${d} already present`); continue; }
  if (!existsSync(src)) { console.error(`no ${d} in main checkout; run ${fetchHint}`); continue; }
  symlinkSync(src, dst, 'junction'); // type is ignored outside Windows
  console.log(`linked ${d}`);
}

function safeLstat(p) {
  try { return lstatSync(p); } catch { return null; }
}
