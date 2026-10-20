#!/usr/bin/env bash
# Writes the golden reference of the old MetaboliteDetector 3.3 for the example runs.
# Run on Linux/WSL with the old binary from Old/MetaboliteDetector-3.3.20200313-Linux_0.deb installed
# (or MD_APP pointing to metabolitedetector_app), g++, nm and python3:
#   docs/reverse-engineering/tools/golden/make-reference.sh
# Output: crates/md-core/tests/data/old-md-3.3.tsv (see golden.py for the format).
# With MD_RESULTS=DIR it also saves DIR/<run>.cmp per run, a known calibration as DIR/harness.calib
# (text) and DIR/harness.calib.bin (binary), and the raw dump DIR/raw.tsv. These come from the old
# program: keep them local (example-data/old-results/ is ignored by git).
# With --reread IN OUT [IN OUT]... it only reads each .cmp, .calib or .calib.bin with the old code
# and writes it again to OUT (checks that the old program accepts files written by md-core).
set -euo pipefail
here=$(cd "$(dirname "$0")" && pwd)
repo=$(cd "$here/../../../.." && pwd)
app=${MD_APP:-/usr/bin/metabolitedetector_app}
data=${MD_DATA:-$repo/example-data/TargetSearchData/inst/gc-ms-data}
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
    want["_ZN4gcms12GCMSDiskScanIifEC1EPKcbbPNS_16ProgressListenerE"] = "DISK_SCAN_CTOR"
    want["_ZNK4gcms12GCMSDiskScanIifE12getScanCountEv"] = "DISK_SCAN_COUNT"
    want["_ZNK4gcms12GCMSDiskScanIifE16getRetentionTimeEj"] = "DISK_SCAN_RT"
    want["_ZN4gcms19DefaultPeakDetectorIifEC1Ev"] = "DEFAULT_DETECTOR_CTOR"
    want["_ZN4gcms12GCMSAnalyzerIifEC1ERKNS_12GCMSDiskScanIifEEPNS_20AbstractPeakDetectorIifEE"] = "ANALYZER_CTOR"
    want["_ZN4gcms12GCMSAnalyzerIifE14detectAllPeaksEbRKfS3_"] = "DETECT_ALL_PEAKS"
    want["_ZN4gcms12GCMSAnalyzerIifE23deconvoluteChromatogramEv"] = "DECONVOLUTE"
    want["_ZNK4gcms13PeakEssentialIfE6getIonEv"] = "PE_ION"
    want["_ZNK4gcms13PeakEssentialIfE12getFirstScanEv"] = "PE_FIRST"
    want["_ZNK4gcms13PeakEssentialIfE11getLastScanEv"] = "PE_LAST"
    want["_ZNK4gcms13PeakEssentialIfE10getMaxScanEv"] = "PE_MAX"
    want["_ZN4gcms13PeakEssentialIfE22getInterpolatedMaxScanEv"] = "PE_INTERP"
    want["_ZNK4gcms13PeakEssentialIfE11getIntegralEv"] = "PE_INTEGRAL"
    want["_ZNK4gcms13PeakEssentialIfE10getQualityEv"] = "PE_QUALITY"
    want["_ZNK4gcms13PeakEssentialIfE12getSharpnessEv"] = "PE_SHARPNESS"
    want["_ZNK4gcms13PeakEssentialIfE16getSignalToNoiseEv"] = "PE_SN"
    want["_ZNK4gcms13PeakEssentialIfE17getModelParameterEv"] = "PE_MODEL"
    want["_ZNK4gcms8SpectrumIifE16getRetentionTimeEv"] = "SPEC_RT"
    want["_ZNK4gcms8SpectrumIifE7getScanEv"] = "SPEC_SCAN"
    want["_ZNK4gcms8SpectrumIifE9getMassesEv"] = "SPEC_MASSES"
    want["_ZNK4gcms8SpectrumIifE9getValuesEv"] = "SPEC_VALUES"
    want["_ZNK4gcms8CompoundIifE17getPeakEssentialsEv"] = "COMPOUND_PEAKS"
    want["_ZN4gcms8CompoundIifE6toDiskERKSt6vectorIPS1_SaIS3_EEPKc"] = "COMPOUND_TO_DISK"
    want["_ZN4gcms8CompoundIifE8fromDiskEPKc"] = "COMPOUND_FROM_DISK"
    want["_ZN4gcms12RICalculatorC1ERKNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEE"] = "RI_CTOR"
    want["_ZN4gcms12RICalculator18addCalibrationPairERKdS2_"] = "RI_ADD_PAIR"
    want["_ZN4gcms12RICalculator6toDiskEPKc"] = "RI_TO_DISK"
    want["_ZN4gcms12RICalculator8fromDiskEPKc"] = "RI_FROM_DISK"
    want["_ZNK4gcms12RICalculator8toStreamERSo"] = "RI_TO_STREAM"
    want["_ZN4gcms12RICalculator10fromStreamERSi"] = "RI_FROM_STREAM"
  }
  ($3 in want) && !seen[$3]++ { print "SYM(" want[$3] ", 0x" $1 ")"; n++ }
  END { if (n != length(want)) { print "missing symbols" > "/dev/stderr"; exit 1 } }
' > "$work/syms.inc"

g++ -std=c++17 -O1 -shared -fPIC -I"$work" -o "$work/old_md_dump.so" "$here/old_md_dump.cpp" -ldl
if [ "${1:-}" = --reread ]; then
  LD_PRELOAD="$work/old_md_dump.so" "$app" "$@"
  exit
fi
LD_PRELOAD="$work/old_md_dump.so" "$app" "$work" "$work/raw.tsv" "$data"/*.cdf > "$work/old.log" 2>&1 || { tail "$work/old.log"; exit 1; }
python3 "$here/golden.py" < "$work/raw.tsv" > "$repo/crates/md-core/tests/data/old-md-3.3.tsv"
# The raw dump has every peak field with 9 digits; the .cmp test compares against it.
if [ -n "${MD_RESULTS:-}" ]; then cp "$work/raw.tsv" "$MD_RESULTS/raw.tsv"; fi
