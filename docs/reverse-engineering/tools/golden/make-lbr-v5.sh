#!/usr/bin/env bash
# Makes version 5 libraries with the old MetaboliteDetector 3.3 and what it reads back from them.
# Run on Linux/WSL with the old binary from Old/MetaboliteDetector-3.3.20200313-Linux_0.deb installed
# (or MD_APP pointing to metabolitedetector_app), g++ and nm:
#   docs/reverse-engineering/tools/golden/make-lbr-v5.sh OUT_DIR [RUN.cdf]
# OUT_DIR gets lib_aa_v5.lbr, fragments_v5.lbr, lib_aa.msl, msl_v5.lbr, run_v5.lbr (with RUN.cdf) and
# old.tsv. Never commit them. Check them with
#   MD_LBR_V5_DIR=OUT_DIR MD_LBR_SAMPLE=lib_aa.lbr cargo test -p md-core --test lbr_sample
set -euo pipefail
here=$(cd "$(dirname "$0")" && pwd)
out=$(mkdir -p "$1" && cd "$1" && pwd)
app=${MD_APP:-/usr/bin/metabolitedetector_app}
lib=${MD_LBR_SAMPLE:-/usr/share/metabolitedetector/lib_aa.lbr}
work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT

# The driver calls into the binary by address; refuse any other build of it.
echo "7883d70cb070c4bd5dc3dc672a914d68d4708cbeecff8c97fe7e5d53bcc64794  $app" | sha256sum -c --quiet

nm "$app" | awk '
  BEGIN {
    want["main"] = "MAIN"
    want["_ZN4gcms15LibraryCompoundIifE8fromDiskEPKc"] = "LC_FROM_DISK"
    want["_ZN4gcms15LibraryCompoundIifE6toDiskERKSt6vectorIPS1_SaIS3_EEPKc"] = "LC_TO_DISK"
    want["_ZN4gcms15LibraryCompoundIifE11addFragmentEmmRKNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEEm"] = "LC_ADD_FRAGMENT"
    want["_ZN4gcms13LibrarySearchIifEC1ERKSt6vectorIPNS_15LibraryCompoundIifEESaIS5_EEb"] = "LS_CTOR_LIB"
    want["_ZN4gcms13LibrarySearchIifEC1ERKSt6vectorIPNS_8CompoundIifEESaIS5_EE"] = "LS_CTOR_COMPOUNDS"
    want["_ZNK4gcms13LibrarySearchIifE9exportMSLEPKcNS1_10ExportTypeE"] = "LS_EXPORT_MSL"
    want["_ZN4gcms13LibrarySearchIifE9importMSLEPKc"] = "LS_IMPORT_MSL"
    want["_ZNK4gcms13LibrarySearchIifE6toDiskEPKc"] = "LS_TO_DISK"
    want["_ZNK4gcms8CompoundIifE17getPeakEssentialsEv"] = "COMPOUND_PEAKS"
    want["_ZN4gcms11CdfImporterC1EPKc"] = "CDF_IMPORTER_CTOR"
    want["_ZN4gcms11CdfImporter10importDataEPKc"] = "CDF_IMPORT_DATA"
    want["_ZN4gcms11CdfImporterD1Ev"] = "CDF_IMPORTER_DTOR"
    want["_ZN4gcms12GCMSDiskScanIifEC1EPKcbbPNS_16ProgressListenerE"] = "DISK_SCAN_CTOR"
    want["_ZN4gcms19DefaultPeakDetectorIifEC1Ev"] = "DEFAULT_DETECTOR_CTOR"
    want["_ZN4gcms12GCMSAnalyzerIifEC1ERKNS_12GCMSDiskScanIifEEPNS_20AbstractPeakDetectorIifEE"] = "ANALYZER_CTOR"
    want["_ZN4gcms12GCMSAnalyzerIifE14detectAllPeaksEbRKfS3_"] = "DETECT_ALL_PEAKS"
    want["_ZN4gcms12GCMSAnalyzerIifE23deconvoluteChromatogramEv"] = "DECONVOLUTE"
  }
  ($3 in want) && !seen[$3]++ { print "SYM(" want[$3] ", 0x" $1 ")"; n++ }
  END { if (n != length(want)) { print "missing symbols" > "/dev/stderr"; exit 1 } }
' > "$work/syms.inc"

g++ -std=c++17 -O1 -shared -fPIC -I"$work" -o "$work/lbr_v5.so" "$here/lbr_v5.cpp" -ldl
# The CDF import writes its .bin/.idx next to OUT_DIR/run; they are removed afterwards.
LD_PRELOAD="$work/lbr_v5.so" "$app" "$out" "$lib" ${2:+"$2"} > "$work/old.log" 2>&1 || { tail "$work/old.log"; exit 1; }
rm -f "$out"/run.bin "$out"/run.idx
ls -l "$out"
