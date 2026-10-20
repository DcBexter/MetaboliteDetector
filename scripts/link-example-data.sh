#!/usr/bin/env bash
# Links the main checkout's example data into this worktree (symlinks), so golden tests do not skip.
set -eu
root=$(git rev-parse --show-toplevel)
main=$(dirname "$(cd "$(git rev-parse --git-common-dir)" && pwd)")
for d in TargetSearchData old-results; do
  dst="$root/example-data/$d"
  src="$main/example-data/$d"
  if [ -e "$dst" ]; then echo "ok: $d already present"; continue; fi
  if [ ! -e "$src" ]; then echo "no $d in main checkout; run example-data/fetch.sh" >&2; continue; fi
  ln -s "$src" "$dst"
  echo "linked $d"
done
