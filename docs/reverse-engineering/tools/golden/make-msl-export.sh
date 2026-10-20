#!/usr/bin/env bash
# Exports a library with the old MetaboliteDetector 3.3 as AMDIS MSL (the three ExportType values).
# Run on Linux/WSL with the old binary installed (or MD_APP pointing to metabolitedetector_app), g++ and nm:
#   docs/reverse-engineering/tools/golden/make-msl-export.sh LIB.lbr|LIB.msl OUT_PREFIX
# A .msl input goes through the old MSL import first (to export libraries with fractional masses).
# Output: OUT_PREFIX.profile.msl, OUT_PREFIX.nominal.msl, OUT_PREFIX.fwhm.msl. Do not commit them.
set -euo pipefail
here=$(cd "$(dirname "$0")" && pwd)
app=${MD_APP:-/usr/bin/metabolitedetector_app}
work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT

echo "7883d70cb070c4bd5dc3dc672a914d68d4708cbeecff8c97fe7e5d53bcc64794  $app" | sha256sum -c --quiet

nm "$app" | awk '
  BEGIN {
    want["main"] = "MAIN"
    want["_ZN4gcms13LibrarySearchIifE8fromDiskEPKc"] = "LIBRARY_SEARCH_FROM_DISK"
    want["_ZN4gcms13LibrarySearchIifE9importMSLEPKc"] = "LIBRARY_SEARCH_IMPORT_MSL"
    want["_ZNK4gcms13LibrarySearchIifE9exportMSLEPKcNS1_10ExportTypeE"] = "LIBRARY_SEARCH_EXPORT_MSL"
  }
  ($3 in want) && !seen[$3]++ { print "SYM(" want[$3] ", 0x" $1 ")"; n++ }
  END { if (n != length(want)) { print "missing symbols" > "/dev/stderr"; exit 1 } }
' > "$work/syms.inc"

g++ -std=c++17 -O1 -shared -fPIC -I"$work" -o "$work/old_msl_export.so" "$here/old_msl_export.cpp" -ldl
LD_PRELOAD="$work/old_msl_export.so" "$app" "$1" "$2"
