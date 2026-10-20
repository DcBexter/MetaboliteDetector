#!/bin/sh
# Runs fixed and random cases through tools::Statistics::oneWayANOVA of the old binary (old_anova.cpp).
# usage: make-anova.sh OUT_DIR [APP]   (in WSL/Linux with g++ and python3)
# Writes OUT_DIR/cases.txt and OUT_DIR/old.tsv (never commit them), then check with
#   MD_ANOVA_DIR=OUT_DIR cargo test -p md-core --test old_anova
set -eu
here=$(cd "$(dirname "$0")" && pwd)
out=$1
app=${2:-$HOME/re/metabolitedetector_app}
echo "7883d70cb070c4bd5dc3dc672a914d68d4708cbeecff8c97fe7e5d53bcc64794  $app" | sha256sum -c -
mkdir -p "$out"
g++ -O1 -std=c++17 -shared -fPIC -o "$out/old_anova.so" "$here/old_anova.cpp" -ldl
python3 "$here/gen_anova.py" > "$out/cases.txt"
# timeout outside the preload, else its main() is replaced too.
timeout 60 env LD_PRELOAD="$out/old_anova.so" "$app" "$out/cases.txt" 2>"$out/err.txt" | grep -E '^(anova|error)	' > "$out/old.tsv"
wc -l "$out/old.tsv"
