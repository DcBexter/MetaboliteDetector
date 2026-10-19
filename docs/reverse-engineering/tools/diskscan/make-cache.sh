#!/usr/bin/env bash
# Converts the example runs into the old MetaboliteDetector 3.3 .bin/.idx cache, for the test
# old_cache_matches_cdf (crates/md-core/tests/example_data.rs). Never commit the output.
# Run on Linux/WSL with the old binary installed (or MD_APP pointing to it), g++ and nm:
#   docs/reverse-engineering/tools/diskscan/make-cache.sh OUT_DIR [CDF...]
#   MD_OLD_CACHE=OUT_DIR cargo test -p md-core --test example_data old_cache
set -euo pipefail
here=$(cd "$(dirname "$0")" && pwd)
repo=$(cd "$here/../../../.." && pwd)
app=${MD_APP:-/usr/bin/metabolitedetector_app}
out=$(mkdir -p "$1" && cd "$1" && pwd)
shift
[[ $# -gt 0 ]] || set -- "$repo"/example-data/TargetSearchData/inst/gc-ms-data/*.cdf
work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT

# The driver calls into the binary by address; refuse any other build of it.
echo "7883d70cb070c4bd5dc3dc672a914d68d4708cbeecff8c97fe7e5d53bcc64794  $app" | sha256sum -c --quiet

nm "$app" | awk '
  BEGIN {
    want["main"] = "MAIN"
    want["_ZN4gcms11CdfImporterC1EPKc"] = "CDF_IMPORTER_CTOR"
    want["_ZN4gcms11CdfImporter10importDataEPKc"] = "CDF_IMPORT_DATA"
    want["_ZN4gcms11CdfImporterD1Ev"] = "CDF_IMPORTER_DTOR"
  }
  ($3 in want) && !seen[$3]++ { print "SYM(" want[$3] ", 0x" $1 ")"; n++ }
  END { if (n != length(want)) { print "missing symbols" > "/dev/stderr"; exit 1 } }
' > "$work/syms.inc"

g++ -std=c++17 -O1 -shared -fPIC -I"$work" -o "$work/old_md_import.so" "$here/old_md_import.cpp" -ldl
LD_PRELOAD="$work/old_md_import.so" "$app" "$out" "$@"
