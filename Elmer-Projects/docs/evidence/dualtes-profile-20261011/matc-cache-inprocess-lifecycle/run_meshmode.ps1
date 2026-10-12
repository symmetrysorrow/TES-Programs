param(
    [Parameter(Mandatory=$true)][string]$ElmerHome,
    [string]$Mpiexec = "C:\Program Files\Microsoft MPI\Bin\mpiexec.exe",
    [string]$MpiBin = "C:\Program Files\Microsoft MPI\Bin",
    [string]$DependencyBin = "G:\Elmer-Migration\elmer\bin"
)
$ErrorActionPreference = "Stop"
$Here = Split-Path -Parent $MyInvocation.MyCommand.Path
$Logs = Join-Path $Here "logs"
$Solver = Join-Path $ElmerHome "bin/ElmerSolver_mpi.exe"
if (-not (Test-Path $Solver)) { throw "ElmerSolver_mpi.exe not found under $ElmerHome" }
if (-not (Test-Path $Mpiexec)) { throw "mpiexec not found: $Mpiexec" }
New-Item -ItemType Directory -Force $Logs | Out-Null
$env:ELMER_HOME = $ElmerHome
$env:OMP_NUM_THREADS = "1"
$env:OMP_DYNAMIC = "FALSE"
$env:PATH = "$(Join-Path $ElmerHome 'bin');$(Join-Path $ElmerHome 'share/elmersolver/lib');$MpiBin;$DependencyBin;G:\Elmer-Migration\msys64\ucrt64\bin;G:\Elmer-Migration\msys64\usr\bin;$env:PATH"

Push-Location $Here
try {
    foreach ($mode in @("on", "off")) {
        $sif = Join-Path $Here "meshmode_cache_$mode.sif"
        $log = Join-Path $Logs "cache_${mode}_raw.log"
        & $Mpiexec -n 1 $Solver $sif 2>&1 | Out-File -LiteralPath $log -Encoding utf8
        $solverExit = $LASTEXITCODE
        if ($solverExit -ne 0) { throw "Cache $mode run failed with exit code $solverExit; see $log" }
    }
    py -3 .\check_lifecycle.py
    if ($LASTEXITCODE -ne 0) { throw "lifecycle evidence checks failed; see pass_fail_skip.json" }
} finally {
    Pop-Location
}
