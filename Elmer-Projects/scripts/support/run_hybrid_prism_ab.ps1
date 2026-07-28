[CmdletBinding()]
param(
    [string]$StandardSolver = 'C:\Program Files\Elmer 26.1-Release\bin\ElmerSolver.exe',
    [string]$CustomSolver = (Join-Path $PSScriptRoot '..\..\..\tools\elmer-hypre\install\bin\ElmerSolver.exe'),
    [string]$StandardRuntimeBin = 'C:\Program Files\Elmer 26.1-Release\bin',
    [string]$CustomRuntimeBin = 'C:\msys64\ucrt64\bin',
    [string]$CaseSif = 'generated\cases\case_tes_steady_hybrid_prism.sif',
    [int]$TimeoutSeconds = 900,
    [switch]$CustomOnly,
    [ValidateSet('Original', 'NoUdfZero')]
    [string]$HeatSourceVariant = 'Original',
    [string]$UdfDll
)

# Copies inputs into a timestamped artifact tree; it never writes the source mesh,
# generated SIF, or either Elmer installation.  The copied UDF is required because
# this steady SIF calls TESTransientHeatSource.
$ErrorActionPreference = 'Stop'
$project = (Resolve-Path (Join-Path $PSScriptRoot '..\..')).Path
if (-not $UdfDll) { $UdfDll = Join-Path $project 'tes_transient_heat_source_t0.dll' }
$stamp = Get-Date -Format 'yyyyMMdd-HHmmss'
$artifactRoot = Join-Path $project "artifacts\hybrid_prism_diagnostics\$stamp"
$sourceSif = Join-Path $project $CaseSif
if (-not (Test-Path -LiteralPath $sourceSif)) { throw "SIF not found: $sourceSif" }
$sifText = Get-Content -LiteralPath $sourceSif -Raw
if ($sifText -notmatch 'Mesh DB "\." "([^"]+)"') { throw 'Could not locate Mesh DB in SIF.' }
$meshName = $Matches[1]
$sourceMesh = Join-Path $project $meshName
if (-not (Test-Path -LiteralPath (Join-Path $sourceMesh 'mesh.header'))) { throw "Mesh not found: $sourceMesh" }

$diag = Join-Path $project 'scripts\support\diagnose_elmer_runtime.py'
function Invoke-IsolatedCase([string]$Label, [string]$Solver, [string]$RuntimeBin) {
    $runDir = Join-Path $artifactRoot $Label
    New-Item -ItemType Directory -Force -Path $runDir | Out-Null
    $oldPath = $env:PATH; $oldElmerHome = $env:ELMER_HOME
    $exitCode = $null; $timedOut = $false; $failure = $null; $solverPath = $Solver; $prefix = $null; $udf = $null
    try {
        $solverPath = (Resolve-Path -LiteralPath $Solver).Path
        Copy-Item -LiteralPath $sourceMesh -Destination (Join-Path $runDir $meshName) -Recurse
        $isolatedSif = Join-Path $runDir 'case.sif'
        Copy-Item -LiteralPath $sourceSif -Destination $isolatedSif
        if ($HeatSourceVariant -eq 'NoUdfZero') {
            $caseText = Get-Content -LiteralPath $isolatedSif -Raw
            $pattern = '(?ms)^  Volumetric Heat Source = Variable Temperature\r?\n    Real Procedure "tes_transient_heat_source_t0" "TESTransientHeatSource"'
            $caseText = [regex]::Replace($caseText, $pattern, '  Volumetric Heat Source = Real 0.0', 1)
            if ($caseText -match 'TESTransientHeatSource') { throw 'NoUdfZero substitution did not remove TESTransientHeatSource.' }
            Set-Content -LiteralPath $isolatedSif -Value $caseText -NoNewline -Encoding utf8
        } elseif ($HeatSourceVariant -eq 'Original') {
            $udf = (Resolve-Path -LiteralPath $UdfDll).Path
            # Procedure names are resolved by this fixed DLL basename.
            Copy-Item -LiteralPath $udf -Destination (Join-Path $runDir 'tes_transient_heat_source_t0.dll')
        }
        $prefix = Split-Path (Split-Path $solverPath -Parent) -Parent
        # The diagnostic must see exactly the loader search order of the child.
        # In particular, custom builds must not silently borrow gfortran DLLs
        # from the standard Elmer installation inherited through PATH.
        $env:ELMER_HOME = $prefix
        $env:PATH = "$(Split-Path $solverPath -Parent);$prefix\share\elmersolver\lib;$RuntimeBin;$oldPath"
        $diagArgs = @($diag, '--solver', $solverPath, '--prefix', $prefix, '--label', $Label)
        if ($HeatSourceVariant -eq 'Original') { $diagArgs += @('--extra-dll', $udf) }
        $diagArgs += @('--output', (Join-Path $runDir 'runtime_manifest.json'))
        & python @diagArgs | Out-Host
        if ($LASTEXITCODE) { throw "runtime diagnostic failed with exit $LASTEXITCODE" }
        $stdout = Join-Path $runDir 'stdout.log'; $stderr = Join-Path $runDir 'stderr.log'
        $process = Start-Process -FilePath $solverPath -ArgumentList 'case.sif' -WorkingDirectory $runDir `
            -RedirectStandardOutput $stdout -RedirectStandardError $stderr -PassThru
        if (-not $process.WaitForExit($TimeoutSeconds * 1000)) {
            Stop-Process -Id $process.Id -Force
            $process.WaitForExit(); $timedOut = $true
        }
        $exitCode = $process.ExitCode
    } catch { $failure = $_ | Out-String } finally {
        $env:PATH = $oldPath
        if ($null -eq $oldElmerHome) { Remove-Item Env:ELMER_HOME -ErrorAction SilentlyContinue } else { $env:ELMER_HOME = $oldElmerHome }
        $udfHash = if ($udf) { (Get-FileHash -LiteralPath $udf -Algorithm SHA256).Hash } else { $null }
        [ordered]@{ label=$Label; solver=$solverPath; prefix=$prefix; heat_source_variant=$HeatSourceVariant; udf_source=$udf; udf_sha256=$udfHash; exit_code=$exitCode; timed_out=$timedOut; failure=$failure; run_dir=$runDir } |
            ConvertTo-Json | Set-Content -LiteralPath (Join-Path $runDir 'run_manifest.json') -Encoding utf8
    }
    [pscustomobject]@{ Label=$Label; ExitCode=$exitCode; TimedOut=$timedOut; Failure=$failure; RunDir=$runDir }
}

New-Item -ItemType Directory -Force -Path $artifactRoot | Out-Null
$outcomes = @()
if (-not $CustomOnly) { $outcomes += Invoke-IsolatedCase 'standard' $StandardSolver $StandardRuntimeBin }
$outcomes += Invoke-IsolatedCase 'custom' $CustomSolver $CustomRuntimeBin
Write-Host "A/B artifacts: $artifactRoot"
if ($outcomes | Where-Object { $_.Failure -or $_.TimedOut -or $_.ExitCode -ne 0 }) { exit 1 }
