[CmdletBinding()]
param(
    [ValidateSet('smoke','full1000')][string]$Mode = 'smoke',
    [string]$MatlabExe = 'matlab'
)
Set-StrictMode -Version 2.0
$ErrorActionPreference = 'Stop'
$Root = Split-Path -Parent $MyInvocation.MyCommand.Path
$Man = Get-Content -LiteralPath (Join-Path $Root 'inputs\required_binary_inputs.json') -Raw | ConvertFrom-Json
foreach ($Item in $Man) {
    $P = Join-Path $Root $Item.destination
    if (-not (Test-Path -LiteralPath $P -PathType Leaf)) { throw "Missing required local input: $P. See docs/DATA_REQUIREMENTS.md." }
    if ((Get-FileHash -LiteralPath $P -Algorithm SHA256).Hash -ne $Item.sha256) { throw "Input hash mismatch: $P" }
}
$EscapedRoot = $Root.Replace('\','/').Replace("'","''")
$Expr = "addpath('$EscapedRoot/reproduction'); run_clean_replay('$EscapedRoot','$Mode');"
$OriginalMatlabPath = $env:MATLABPATH
try {
    $env:MATLABPATH = ''
    & $MatlabExe -wait -sd $Root -batch $Expr
    $RC = $LASTEXITCODE
    if ($RC -ne 0) { throw "MATLAB replay failed (exit=$RC). Inspect validation/replay_status.json." }
} finally {
    $env:MATLABPATH = $OriginalMatlabPath
}
