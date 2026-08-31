param(
    [ValidateSet("original", "hybrid")]
    [string]$TimeGrid = "original",
    # One RTX 3060 Ti is normally fastest with one AMGX rank.  Override with
    # -MpiProcs 4 only when deliberately comparing against the four-rank run.
    [int]$MpiProcs = 1,
    [int]$SmokeSteps = 0,
    [string]$Distro = "Ubuntu",
    [string]$AmgxConfig = "config\amgx\tes_fgmres_aggregation_l1_5e-8.json",
    [ValidateSet("default", "slave", "master", "slave-transpose", "master-transpose", "dual-lagrange")]
    [string]$AmgxConstraintMode = "default",
    [switch]$ReuseSteady,
    [switch]$ReuseKnownSteady,
    [switch]$ForceDeps,
    [switch]$DryRun
)

$ErrorActionPreference = "Stop"

if ($MpiProcs -lt 1) {
    throw "MpiProcs must be at least 1."
}
if ($SmokeSteps -lt 0) {
    throw "SmokeSteps must be zero (disabled) or a positive integer."
}
if ($ReuseSteady -and $ReuseKnownSteady) {
    throw "ReuseSteady and ReuseKnownSteady are mutually exclusive."
}
if ($ForceDeps -and $ReuseKnownSteady) {
    throw "ForceDeps cannot be combined with a preexisting known steady result."
}

$repo = (Resolve-Path (Join-Path $PSScriptRoot "..")).Path
$tools = Join-Path (Split-Path $repo -Parent) "tools"
$solverWin = Join-Path $tools "elmer-gpu-wsl\bin\ElmerSolver_mpi"
$amgxWin = Join-Path $repo $AmgxConfig

foreach ($path in @($solverWin, $amgxWin)) {
    if (-not (Test-Path -LiteralPath $path -PathType Leaf)) {
        throw "Required GPU file is missing: $path"
    }
}

function Convert-ToWslPath([string]$WindowsPath) {
    # wsl.exe removes backslashes while parsing a direct wslpath argument in
    # PowerShell.  These project paths are on a normal local drive, so the
    # explicit /mnt/<drive> conversion is both reliable and transparent.
    if ($WindowsPath -notmatch '^[A-Za-z]:\\') {
        throw "Expected a Windows drive path: $WindowsPath"
    }
    $drive = $WindowsPath.Substring(0, 1).ToLowerInvariant()
    $rest = $WindowsPath.Substring(2) -replace '\\', '/'
    return "/mnt/$drive$rest"
}

function Bash-Quote([string]$Value) {
    # The current repository is on a path without spaces.  Keep these command
    # fragments unquoted so PowerShell does not reinterpret nested bash quotes.
    return $Value
}

$repoWsl = Convert-ToWslPath $repo
$toolsWsl = Convert-ToWslPath $tools
$solverWsl = "$toolsWsl/elmer-gpu-wsl/bin/ElmerSolver_mpi"
$amgxWsl = Convert-ToWslPath $amgxWin
$amgxLibWsl = "$toolsWsl/amgx-gpu-install-mpi/lib"
$elmerWsl = "$toolsWsl/elmer-gpu-wsl"
$fmodulesWsl = "/home/symme/elmer-gpu-build/fmodules"
$elmerBuildLibWsl = "/home/symme/elmer-gpu-build/fem/src"
$parallelUdfWsl = "$repoWsl/tes_parallel_circuit.so"
$transientUdfWsl = "$repoWsl/tes_transient_heat_source_t0.so"

$buildUdfs = @(
    "test -f $(Bash-Quote "$fmodulesWsl/defutils.mod")",
    "gfortran -O2 -fPIC -shared -I$(Bash-Quote $fmodulesWsl) $(Bash-Quote "$repoWsl/tes_parallel_circuit.f90") -L$(Bash-Quote $elmerBuildLibWsl) -L$(Bash-Quote $amgxLibWsl) -Wl,-rpath,$(Bash-Quote "$elmerWsl/lib/elmersolver") -Wl,-rpath,$(Bash-Quote $amgxLibWsl) -lelmersolver -lamgxsh -o $(Bash-Quote $parallelUdfWsl)",
    "gfortran -O2 -fPIC -shared -I$(Bash-Quote $fmodulesWsl) $(Bash-Quote "$repoWsl/tes_transient_heat_source.f90") -L$(Bash-Quote $elmerBuildLibWsl) -L$(Bash-Quote $amgxLibWsl) -Wl,-rpath,$(Bash-Quote "$elmerWsl/lib/elmersolver") -Wl,-rpath,$(Bash-Quote $amgxLibWsl) -lelmersolver -lamgxsh -o $(Bash-Quote $transientUdfWsl)"
)

$arguments = @(
    "python3 scripts/prep/run_singlepixel_prod_v2_original_timegrid.py",
    "--mpi-procs $MpiProcs",
    "--time-grid $TimeGrid",
    "--linear-system mumps",
    "--elmer-solver $(Bash-Quote $solverWsl)",
    "--runtime-bin ''",
    "--amgx-config $(Bash-Quote $amgxWsl)"
    "--amgx-constraint-mode $AmgxConstraintMode"
)
if ($ReuseSteady) {
    $arguments += "--reuse-steady"
}
if ($ReuseKnownSteady) {
    $arguments += "--reuse-known-serial-steady"
}
if ($ForceDeps) {
    $arguments += "--force-deps"
}
if ($SmokeSteps -gt 0) {
    $arguments += "--smoke-steps $SmokeSteps"
}
if ($DryRun) {
    $arguments += "--dry-run"
}

$bashScript = @(
    "set -e",
    "export ELMER_HOME=$(Bash-Quote $elmerWsl)",
    "export PATH=$(Bash-Quote "$elmerWsl/bin"):/usr/lib/wsl/lib:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin",
    # CUDA must load WSL's paravirtualized driver before any native Linux
    # libcuda package that may also be installed in the distro.
    "export LD_LIBRARY_PATH=/usr/lib/wsl/lib:$(Bash-Quote $repoWsl):$(Bash-Quote $amgxLibWsl):$(Bash-Quote "$elmerWsl/lib/elmersolver"):/usr/lib/x86_64-linux-gnu:/usr/lib/x86_64-linux-gnu/openmpi/lib:`$LD_LIBRARY_PATH",
    "cd $(Bash-Quote $repoWsl)",
    ($buildUdfs -join "; `n"),
    ($arguments -join " ")
) -join "; `n"

Write-Host "Starting GPU Elmer through WSL ($Distro):"
Write-Host ($arguments -join " ")
& wsl.exe -d $Distro -- bash -lc $bashScript
exit $LASTEXITCODE
