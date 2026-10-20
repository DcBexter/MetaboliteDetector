#!/usr/bin/env bash
# Runs the old binary's centroiding and binning on a cases file, for the test matches_old_program
# (crates/md-core/tests/centroid.rs). Never commit the output. Run on Linux/WSL with the old binary
# installed (or MD_APP pointing to it), g++ and nm:
#   MD_SPEC_CASES=DIR cargo test -p md-core --test centroid write_cases
#   docs/reverse-engineering/tools/spectrum/run.sh DIR      # reads DIR/cases.txt, writes DIR/old.txt
#   MD_SPEC_OLD=DIR cargo test -p md-core --test centroid matches_old_program
set -euo pipefail
here=$(cd "$(dirname "$0")" && pwd)
app=${MD_APP:-/usr/bin/metabolitedetector_app}
dir=$(cd "$1" && pwd)
work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT

# The driver calls into the binary by address; refuse any other build of it.
echo "7883d70cb070c4bd5dc3dc672a914d68d4708cbeecff8c97fe7e5d53bcc64794  $app" | sha256sum -c --quiet

nm "$app" | awk '
  BEGIN {
    want["main"] = "MAIN"
    want["_ZN4gcms8SpectrumIifEC1EPSt6vectorIiSaIiEEPS2_IfSaIfEE"] = "SPECTRUM_CTOR"
    want["_ZN4gcms17SpectrumFunctions19getCentroidSpectrumIifEEPNS_16CentroidSpectrumIT_T0_EERNS_8SpectrumIS3_S4_EERKdSB_SB_"] = "CENTROID"
    want["_ZN4gcms8SpectrumIifE17getBinnedSpectrumERKSt6vectorIiSaIiEERKS2_IfSaIfEEf"] = "BINNED"
  }
  ($3 in want) && !seen[$3]++ { print "SYM(" want[$3] ", 0x" $1 ")"; n++ }
  END { if (n != length(want)) { print "missing symbols" > "/dev/stderr"; exit 1 } }
' > "$work/syms.inc"

g++ -std=c++17 -O1 -shared -fPIC -I"$work" -o "$work/old_spectrum.so" "$here/old_spectrum.cpp" -ldl
LD_PRELOAD="$work/old_spectrum.so" "$app" "$dir/cases.txt" "$dir/old.txt"
