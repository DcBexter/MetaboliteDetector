#!/usr/bin/env bash
# Black-box checks of the NIST text import (SQLLibrary::importNIST) and the library search
# (getPotentialHits, getLibraryHits) of the old MetaboliteDetector 3.3 (Linux). Run on Linux/WSL with
# the old binary from Old/MetaboliteDetector-3.3.20200313-Linux_0.deb installed (or MD_APP set), g++, nm
# and python3 (sqlite3 module):
#   docs/reverse-engineering/tools/nist/make-nist-dump.sh
# All input is synthetic (invented names and peaks, in crates/md-core/tests/data/synthetic-nist.sdf); no NIST data is used or downloaded.
# The script prints the parser fields, the import report, the imported tables and the search results.
# Everything is written to a temp directory that is removed on exit.
set -euo pipefail
here=$(cd "$(dirname "$0")" && pwd)
app=${MD_APP:-/usr/bin/metabolitedetector_app}
work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT

# The driver calls into the binary by address; refuse any other build of it.
echo "7883d70cb070c4bd5dc3dc672a914d68d4708cbeecff8c97fe7e5d53bcc64794  $app" | sha256sum -c --quiet

nm "$app" | awk '
  BEGIN {
    want["main"] = "MAIN"
    want["_ZN4gcms10NISTParserC1EPKc"] = "NISTPARSER_CTOR"
    want["_ZN4gcms10NISTParserD1Ev"] = "NISTPARSER_DTOR"
    want["_ZN4gcms10NISTParser10nextRecordEv"] = "NEXT_RECORD"
    want["_ZN4gcms10NISTParser9nextFieldEv"] = "NEXT_FIELD"
    want["_ZN4gcms10NISTParser12getFieldNameB5cxx11Ev"] = "GET_FIELD_NAME"
    want["_ZN4gcms10NISTParser13getFieldValueB5cxx11Ev"] = "GET_FIELD_VALUE"
    want["_ZN4gcms10SQLLibraryC1EPKc"] = "SQL_CTOR"
    want["_ZN4gcms10SQLLibrary10importNISTEPKc"] = "SQL_IMPORT_NIST"
    want["_ZNK4gcms10SQLLibrary16getPotentialHitsB5cxx11ERKNS_8SpectrumIifEEd"] = "SQL_POTENTIAL"
    want["_ZNK4gcms10SQLLibrary14getLibraryHitsERKNS_8CompoundIifEEf"] = "SQL_LIB_HITS"
    want["_ZN4gcms8SpectrumIifEC1EPSt6vectorIiSaIiEEPS2_IfSaIfEE"] = "SPEC_CTOR"
    want["_ZN4gcms8CompoundIifEC1EdPSt6vectorIiSaIiEEPS2_IfSaIfEEf"] = "COMP_CTOR"
    want["_ZNK4gcms8SpectrumIifE20getSpectrumFragmentsERNSt7__cxx114listINS3_ISt4pairIfiESaIS5_EEESaIS7_EEEdi"] = "SPEC_FRAGMENTS"
    want["_ZN4gcms8SpectrumIifE17setRetentionIndexEd"] = "SPEC_SET_RI"
    want["_ZNK4gcms8CompoundIifE5getIDEv"] = "COMP_ID"
    want["_ZNK4gcms8CompoundIifE7getNameB5cxx11Ev"] = "COMP_NAME"
  }
  ($3 in want) && !seen[$3]++ { print "SYM(" want[$3] ", 0x" $1 ")"; n++ }
  END { if (n != length(want)) { print "missing symbols" > "/dev/stderr"; exit 1 } }
' > "$work/syms.inc"

g++ -std=c++17 -O1 -shared -fPIC -I"$work" -o "$work/nist_dump.so" "$here/nist_dump.cpp" -ldl

# Synthetic NIST-style text (invented names and peaks), shared with the md-core tests. Peaks are
# "mz intensity" lines; the old importer keeps only peaks with mz > 69 and intensity > 0.
nist=$here/../../../../crates/md-core/tests/data/synthetic-nist.sdf

# Adjacent masses (differ by exactly 1.0) form the fragment groups of the binary (Observed, see notes);
# the synthetic spectra include them. Query spectrum: a line "RT <rt>", then "mz intensity" peaks identical to Synthetic alpha.
cat > "$work/query.txt" <<'EOF'
RT 10.5
RI 1455
71 100
72 12
77 40
78 9
91 80
92 15
105 60
106 7
119 50
133 40
147 30
161 20
EOF

# The same query with strong falling pairs below m/z 70, which the library side never stores: shows
# whether the query side drops them before choosing its 5 strongest fragments.
{ head -2 "$work/query.txt"; printf '41 300\n42 30\n43 290\n44 29\n55 280\n56 28\n57 270\n58 27\n'
  tail -n +3 "$work/query.txt"; } > "$work/query-low.txt"

echo "== parse (NISTParser fields) =="
LD_PRELOAD="$work/nist_dump.so" "$app" parse "$nist"

echo "== library (importNIST, getPotentialHits, getLibraryHits) =="
LD_PRELOAD="$work/nist_dump.so" "$app" library "$work/lib.sqlite" "$nist" "$work/query.txt"

echo "== library, query with strong ions below m/z 70 =="
LD_PRELOAD="$work/nist_dump.so" "$app" library "$work/lib-low.sqlite" "$nist" "$work/query-low.txt"

echo "== sqlite tables written by importNIST =="
python3 -I - "$work/lib.sqlite" <<'PY'
import sqlite3, sys
con = sqlite3.connect(sys.argv[1])
for (sql,) in con.execute("SELECT sql FROM sqlite_master WHERE sql IS NOT NULL ORDER BY name"):
    print("schema\t" + " ".join(sql.split()))
for (table,) in con.execute("SELECT name FROM sqlite_master WHERE type='table' ORDER BY name"):
    cur = con.execute('SELECT * FROM "%s"' % table)
    cols = [d[0] for d in cur.description]
    rows = cur.fetchall()
    print("table\t%s\trows\t%d\tcolumns\t%s" % (table, len(rows), ",".join(cols)))
    for row in rows[:40]:
        print("row\t%s\t%s" % (table, "\t".join(repr(v) for v in row)))
PY
