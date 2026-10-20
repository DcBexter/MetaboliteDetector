#!/usr/bin/env bash
# Imports a FastFlight run with the old MetaboliteDetector 3.3 importer and writes its .bin/.idx cache,
# for the test matches_old_importer (crates/md-core/tests/fastflight.rs). Never commit the output.
# Run on Linux/WSL with the old binary installed (or MD_APP pointing to it), g++ and nm:
#   MD_FF_SYNTH=DIR cargo test -p md-core --test fastflight write_synthetic
#   docs/reverse-engineering/tools/fastflight/make-cache.sh DIR
#   MD_FF_OLD=DIR cargo test -p md-core --test fastflight matches_old_importer
# DIR holds MsData.FFT, MsData.FFC and one .mzcalib (calib.mzcalib); the cache goes to DIR/old/synthetic.*
set -euo pipefail
here=$(cd "$(dirname "$0")" && pwd)
app=${MD_APP:-/usr/bin/metabolitedetector_app}
dir=$(cd "$1" && pwd)
mkdir -p "$dir/old"
work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT

# The driver calls into the binary by address; refuse any other build of it.
echo "7883d70cb070c4bd5dc3dc672a914d68d4708cbeecff8c97fe7e5d53bcc64794  $app" | sha256sum -c --quiet

nm "$app" | awk '
  BEGIN {
    want["main"] = "MAIN"
    want["_ZN4gcms18FastFlightImporterC1EPKc"] = "FF_IMPORTER_CTOR"
    want["_ZN4gcms18FastFlightImporter10importDataEPKcbb"] = "FF_IMPORT_DATA"
    want["_ZN4gcms18FastFlightImporter23setCalibrationParameterERKSt6vectorIdSaIdEE"] = "FF_SET_CALIBRATION"
    want["_ZN4gcms18FastFlightImporterD1Ev"] = "FF_IMPORTER_DTOR"
  }
  ($3 in want) && !seen[$3]++ { print "SYM(" want[$3] ", 0x" $1 ")"; n++ }
  END { if (n != length(want)) { print "missing symbols" > "/dev/stderr"; exit 1 } }
' > "$work/syms.inc"

g++ -std=c++17 -O1 -shared -fPIC -I"$work" -o "$work/old_ff_import.so" "$here/old_ff_import.cpp" -ldl
LD_PRELOAD="$work/old_ff_import.so" "$app" "$dir" "$dir/calib.mzcalib" "$dir/old/synthetic"
