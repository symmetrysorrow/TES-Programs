param(
    [Parameter(Mandatory=$true)][string]$ElmerSource,
    [Parameter(Mandatory=$true)][string]$BuildDir,
    [Parameter(Mandatory=$true)][string]$InstallDir,
    [string]$FortranCompiler = "gfortran",
    [string]$CCompiler = "gcc",
    [string]$CxxCompiler = "g++",
    [int]$Parallel = 6
)
$ErrorActionPreference = "Stop"
$Here = Split-Path -Parent $MyInvocation.MyCommand.Path
$Patch = Join-Path $Here "04_matc_cache_meshmode_diagnostics.patch"
if (-not (Test-Path (Join-Path $ElmerSource "fem/src/ElmerSolver.F90"))) { throw "ElmerSource is not an Elmer source tree: $ElmerSource" }
if (-not (Test-Path (Join-Path $ElmerSource ".git"))) { throw "ElmerSource must be a git checkout with patches 01-03 already applied" }

$env:PATH = "G:\Elmer-Migration\msys64\ucrt64\bin;G:\Elmer-Migration\msys64\usr\bin;" + $env:PATH
git -C $ElmerSource apply --check $Patch
if ($LASTEXITCODE -ne 0) { throw "patch 04 does not apply cleanly; expected an unmodified 01-03 source tree" }
git -C $ElmerSource apply $Patch
if ($LASTEXITCODE -ne 0) { throw "failed to apply patch 04" }

cmake -S $ElmerSource -B $BuildDir -G Ninja `
    -DCMAKE_BUILD_TYPE=Release `
    "-DCMAKE_C_COMPILER=$CCompiler" `
    "-DCMAKE_CXX_COMPILER=$CxxCompiler" `
    "-DCMAKE_Fortran_COMPILER=$FortranCompiler" `
    -DCMAKE_Fortran_FLAGS=-DPHASE24_MATC_LIFECYCLE_DIAGNOSTICS `
    "-DCMAKE_INSTALL_PREFIX=$InstallDir" `
    -DELMER_INSTALL_LIB_DIR=share/elmersolver/lib `
    -DWITH_MPI=ON -DWITH_OpenMP=ON -DWITH_Hypre=ON -DWITH_Mumps=OFF `
    -DWITH_ELMERGUI=OFF -DWITH_ELMERPOST=OFF -DWITH_UMFPACK=ON -DWITH_LUA=ON
if ($LASTEXITCODE -ne 0) { throw "CMake configure failed" }
cmake --build $BuildDir --parallel $Parallel
if ($LASTEXITCODE -ne 0) { throw "diagnostic Elmer build failed" }
cmake --install $BuildDir
if ($LASTEXITCODE -ne 0) { throw "diagnostic Elmer install failed" }
Write-Host "Diagnostic solver installed to $InstallDir"
