# Downloads the TargetSearchData example set (Bioconductor, GPL >= 2) into example-data\TargetSearchData.
$ErrorActionPreference = 'Stop'
Set-Location $PSScriptRoot
$url = 'https://bioconductor.org/packages/3.23/data/experiment/src/contrib/TargetSearchData_1.50.0.tar.gz'
$sha256 = '4E42D877FF1EFAFB55183E7BA561D7901672567D28E7D349F9519437FA1C7AB5'
Invoke-WebRequest -Uri $url -OutFile TargetSearchData.tar.gz
if ((Get-FileHash TargetSearchData.tar.gz -Algorithm SHA256).Hash -ne $sha256) {
    Remove-Item TargetSearchData.tar.gz
    throw 'Checksum mismatch, download removed.'
}
New-Item -ItemType Directory -Force TargetSearchData | Out-Null
# tar.exe ships with Windows 10 1803 and later.
tar.exe -xzf TargetSearchData.tar.gz -C TargetSearchData --strip-components=1 `
    TargetSearchData/DESCRIPTION TargetSearchData/README.md TargetSearchData/inst/gc-ms-data
if ($LASTEXITCODE -ne 0) { throw 'Extraction failed.' }
Remove-Item TargetSearchData.tar.gz
Write-Host "Done: $PSScriptRoot\TargetSearchData\inst\gc-ms-data"
