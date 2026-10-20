#!/usr/bin/env bash
# Runs the GSL black-box oracle (gsl_oracle.cpp) inside the old MetaboliteDetector 3.3 binary:
#   docs/reverse-engineering/tools/gsl-oracle/run.sh < commands.txt > outputs.txt
# Linux/WSL with the old .deb installed (or MD_APP pointing to metabolitedetector_app), g++ and nm.
# MD_ORACLE_BUILD=DIR keeps the built driver in DIR between runs. The outputs come from the old
# program: keep them local, do not commit them.
set -euo pipefail
here=$(cd "$(dirname "$0")" && pwd)
app=${MD_APP:-/usr/bin/metabolitedetector_app}
if [[ -n "${MD_ORACLE_BUILD:-}" ]]; then
  work=$MD_ORACLE_BUILD
  mkdir -p "$work"
else
  work=$(mktemp -d)
  trap 'rm -rf "$work"' EXIT
fi

if [[ ! -f "$work/gsl_oracle.so" || "$here/gsl_oracle.cpp" -nt "$work/gsl_oracle.so" ]]; then
  # The oracle calls into the binary by address; refuse any other build of it.
  echo "7883d70cb070c4bd5dc3dc672a914d68d4708cbeecff8c97fe7e5d53bcc64794  $app" | sha256sum -c --quiet
  nm "$app" | awk '
    BEGIN {
      want["main"] = "MAIN"
      want["_ZN5tools10Statistics20polynomialRegressionERSt6vectorIdSaIdEES4_j"] = "POLY_REGRESSION"
      want["gsl_set_error_handler_off"] = "ERROR_HANDLER_OFF"
      want["gsl_matrix_alloc"] = "MATRIX_ALLOC"
      want["gsl_matrix_free"] = "MATRIX_FREE"
      want["gsl_matrix_get"] = "MATRIX_GET"
      want["gsl_matrix_set"] = "MATRIX_SET"
      want["gsl_vector_alloc"] = "VECTOR_ALLOC"
      want["gsl_vector_free"] = "VECTOR_FREE"
      want["gsl_vector_get"] = "VECTOR_GET"
      want["gsl_vector_set"] = "VECTOR_SET"
      want["gsl_multifit_linear_alloc"] = "MFL_ALLOC"
      want["gsl_multifit_linear_free"] = "MFL_FREE"
      want["gsl_multifit_linear"] = "MFL"
      want["gsl_linalg_balance_columns"] = "BALANCE_COLUMNS"
      want["gsl_linalg_SV_decomp"] = "SV_DECOMP"
      want["gsl_linalg_SV_decomp_mod"] = "SV_DECOMP_MOD"
      want["gsl_linalg_SV_solve"] = "SV_SOLVE"
      want["gsl_linalg_bidiag_decomp"] = "BIDIAG_DECOMP"
      want["gsl_linalg_householder_transform"] = "HOUSEHOLDER_TRANSFORM"
      want["gsl_linalg_householder_hm"] = "HOUSEHOLDER_HM"
      want["gsl_linalg_householder_mh"] = "HOUSEHOLDER_MH"
      want["cblas_dnrm2"] = "CBLAS_DNRM2"
      want["cblas_ddot"] = "CBLAS_DDOT"
    }
    ($3 in want) && !seen[$3]++ { print "SYM(" want[$3] ", 0x" $1 ")"; n++ }
    $3 == "xhypot" { print "SYM(XHYPOT" nx++ ", 0x" $1 ")" }
    END { if (n != length(want)) { print "missing symbols" > "/dev/stderr"; exit 1 } }
  ' > "$work/syms.inc"
  g++ -std=c++17 -O1 -shared -fPIC -I"$work" -o "$work/gsl_oracle.so" "$here/gsl_oracle.cpp" -ldl
fi
LD_PRELOAD="$work/gsl_oracle.so" "$app" 2>/dev/null
