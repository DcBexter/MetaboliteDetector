#!/bin/sh
# Downloads the TargetSearchData example set (Bioconductor, GPL >= 2) into example-data/TargetSearchData.
set -eu
cd "$(dirname "$0")"
URL=https://bioconductor.org/packages/3.23/data/experiment/src/contrib/TargetSearchData_1.50.0.tar.gz
SHA256=4e42d877ff1efafb55183e7ba561d7901672567d28e7d349f9519437fa1c7ab5
curl -fL -o TargetSearchData.tar.gz "$URL"
echo "$SHA256  TargetSearchData.tar.gz" | sha256sum -c -
mkdir -p TargetSearchData
tar -xzf TargetSearchData.tar.gz -C TargetSearchData --strip-components=1 \
    TargetSearchData/DESCRIPTION TargetSearchData/README.md TargetSearchData/inst/gc-ms-data
rm TargetSearchData.tar.gz
echo "Done: $(pwd)/TargetSearchData/inst/gc-ms-data"
