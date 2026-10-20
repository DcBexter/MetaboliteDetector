#!/bin/sh
# Runs synthetic MID cases through gcms::MIDCompoundSeries of the old binary (old_midseries.cpp),
# then its group statistics and result window exports (old_midstats.cpp, needs Qt offscreen).
# usage: make-midseries.sh OUT_DIR [APP] [LIBRARY.lbr]   (in WSL/Linux with g++ and python3)
# Writes OUT_DIR/cases.txt, old.tsv and old-stats.txt (never commit them), then check with
#   MD_MID_SERIES_DIR=OUT_DIR MD_LBR_SAMPLE=LIBRARY.lbr cargo test -p md-core --test mid_series
set -eu
here=$(cd "$(dirname "$0")" && pwd)
out=$1
app=${2:-$HOME/re/metabolitedetector_app}
lbr=${3:-$HOME/re/deb/root/usr/share/metabolitedetector/lib_aa.lbr}
echo "7883d70cb070c4bd5dc3dc672a914d68d4708cbeecff8c97fe7e5d53bcc64794  $app" | sha256sum -c -
mkdir -p "$out"
for t in old_midseries old_midstats; do
    g++ -O1 -std=c++17 -shared -fPIC -o "$out/$t.so" "$here/$t.cpp" -ldl
done
python3 "$here/gen_midseries.py" > "$out/cases.txt"
# timeout outside the preload, else its main() is replaced too.
timeout 60 env LD_PRELOAD="$out/old_midseries.so" "$app" "$here/old-isotopes.tsv" "$lbr" "$out/cases.txt" \
    2>"$out/err.txt" | grep -E '^(norm|mid|error)	' > "$out/old.tsv"
wc -l "$out/old.tsv"
QT_QPA_PLATFORM=offscreen timeout 60 env LD_PRELOAD="$out/old_midstats.so" "$app" "$here/old-isotopes.tsv" "$lbr" \
    "$out/cases.txt" 2>"$out/err-stats.txt" > "$out/old-stats.txt"
grep -c '^stat' "$out/old-stats.txt"
