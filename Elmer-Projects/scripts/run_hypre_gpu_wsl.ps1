[CmdletBinding()]
param(
    [ValidateSet("cuda", "hip")]
    [string]$Backend = "cuda",
    [string]$GpuArchitecture = "86",
    [string]$HypreTag = "v3.0.0",
    [string]$Case = "case_p19_hypre_flexgmres_boomeramg_gpu_time5us_smoke_7step",
    [int]$MpiProcs = 1,
    [switch]$Build,
    [switch]$DryRun
)

$ErrorActionPreference = "Stop"
if ($MpiProcs -lt 1) { throw "MpiProcs must be at least one" }

$repo = (Resolve-Path (Join-Path $PSScriptRoot "..")).Path
$tools = Join-Path (Split-Path $repo -Parent) "tools"
$tagSuffix = if ($HypreTag -eq "v3.0.0") { "" } else {
    "-" + (($HypreTag -replace '[^A-Za-z0-9]+', '-').Trim('-'))
}
$prefix = Join-Path $tools "elmer-hypre-$Backend$tagSuffix-wsl"
$solver = Join-Path $prefix "bin\ElmerSolver_mpi"
if ($Build) {
    & (Join-Path $PSScriptRoot "support\build_elmer_hypre_gpu_wsl.ps1") `
        -Backend $Backend -GpuArchitecture $GpuArchitecture -HypreTag $HypreTag
    if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }
}
if (-not (Test-Path -LiteralPath $solver -PathType Leaf)) {
    throw "Missing $solver. Re-run with -Build."
}

python scripts\prep\prepare_hypre_gpu_phase19.py
if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }

function To-WslPath([string]$Path) {
    $drive = $Path.Substring(0, 1).ToLowerInvariant()
    return "/mnt/$drive" + (($Path.Substring(2) -replace '\\', '/'))
}
$repoWsl = To-WslPath $repo
$toolsWsl = To-WslPath $tools
$prefixWsl = To-WslPath $prefix
$solverWsl = "$prefixWsl/bin/ElmerSolver_mpi"
$hypreWsl = "$toolsWsl/hypre-$Backend$tagSuffix-install"
$amgxWsl = "$toolsWsl/amgx-gpu-install-mpi/lib"
$fmodulesWsl = "$prefixWsl/share/elmersolver/include"
$projectWsl = "$repoWsl/elmer_project_hypre_gpu_phase19.json"
$udfCircuit = "$repoWsl/tes_parallel_circuit.so"
$udfPulse = "$repoWsl/tes_transient_heat_source_t0.so"

$runOptions = if ($DryRun) { "--dry-run" } else { "" }
$deviceEnv = if ($Backend -eq "cuda") { "export CUDA_VISIBLE_DEVICES=0" } else { "export HIP_VISIBLE_DEVICES=0" }
$bash = @"
set -euo pipefail
$deviceEnv
export ELMER_HOME='$prefixWsl'
export PATH='$prefixWsl/bin':/usr/lib/wsl/lib:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin
export LD_LIBRARY_PATH='/usr/lib/wsl/lib:${hypreWsl}/lib:${amgxWsl}:${prefixWsl}/lib/elmersolver:${repoWsl}'
cd '$repoWsl'
gfortran -O2 -fPIC -shared -I'$fmodulesWsl' tes_parallel_circuit.f90 -L'$prefixWsl/lib/elmersolver' -L'$amgxWsl' -Wl,-rpath,'$prefixWsl/lib/elmersolver' -Wl,-rpath,'$amgxWsl' -lelmersolver -lamgxsh -o '$udfCircuit'
gfortran -O2 -fPIC -shared -I'$fmodulesWsl' tes_transient_heat_source.f90 -L'$prefixWsl/lib/elmersolver' -L'$amgxWsl' -Wl,-rpath,'$prefixWsl/lib/elmersolver' -Wl,-rpath,'$amgxWsl' -lelmersolver -lamgxsh -o '$udfPulse'
python3 run.py '$Case' --project '$projectWsl' --mpi-procs $MpiProcs --elmer-solver '$solverWsl' --runtime-bin '' $runOptions
"@
Write-Host "Running $Case with HYPRE $Backend ($MpiProcs MPI rank(s))."
& wsl.exe -d Ubuntu -- bash -lc $bash
exit $LASTEXITCODE
