# Links the main checkout's example data into this worktree (directory junctions), so golden tests do not skip.
$ErrorActionPreference = 'Stop'
$root = (git rev-parse --show-toplevel).Trim()
$common = (Resolve-Path (git rev-parse --git-common-dir)).Path
$main = Split-Path $common -Parent
foreach ($d in 'TargetSearchData', 'old-results') {
    $dst = Join-Path $root "example-data/$d"
    $src = Join-Path $main "example-data/$d"
    if (Test-Path $dst) { Write-Host "ok: $d already present"; continue }
    if (-not (Test-Path $src)) { Write-Warning "no $d in main checkout; run example-data/fetch.ps1"; continue }
    New-Item -ItemType Junction -Path $dst -Target $src | Out-Null
    Write-Host "linked $d"
}
