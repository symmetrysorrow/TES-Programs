[CmdletBinding()]
param(
    [string]$Source = (Join-Path $PSScriptRoot '..\..\..\tools\elmer-hypre\src'),
    [string]$Build = (Join-Path $PSScriptRoot '..\..\..\tools\elmer-hypre\build-phase1-gfortran16'),
    [string]$Install = (Join-Path $PSScriptRoot '..\..\..\tools\elmer-hypre\install-phase1-gfortran16'),
    [switch]$DebugBuild
)

# New directory names are intentional. Refuse an existing build/install rather
# than mixing files into a prior runtime. If configure fails, rerun with distinct
# -Build/-Install paths (for example ...-retry1); do not reuse a partial tree.
$ErrorActionPreference = 'Stop'
$sourcePath = (Resolve-Path -LiteralPath $Source).Path
foreach ($path in @($Build, $Install)) { if (Test-Path -LiteralPath $path) { throw "Refusing existing path: $path" } }
$bin = 'C:\msys64\ucrt64\bin'
$cmake = Join-Path $bin 'cmake.exe'; $ninja = Join-Path $bin 'ninja.exe'
foreach ($tool in @($cmake, $ninja, (Join-Path $bin 'gcc.exe'), (Join-Path $bin 'gfortran.exe'), (Join-Path $bin 'mpicc.exe'), (Join-Path $bin 'mpif90.exe'))) { if (-not (Test-Path $tool)) { throw "Tool missing: $tool" } }
$buildType = if ($DebugBuild) { 'Debug' } else { 'Release' }
$fortranFlags = if ($DebugBuild) { '-O0 -g -fcheck=all -fbacktrace' } else { '' }
$openBlas = 'C:\Program Files\Elmer 26.1-Release\bin\libopenblas.dll.a'
$include = Join-Path $bin '..\include'; $lib = Join-Path $bin '..\lib'
$configureArgs = @(
    '-S', $sourcePath, '-B', $Build, '-G', 'Ninja',
    "-DCMAKE_MAKE_PROGRAM=$ninja", "-DCMAKE_BUILD_TYPE=$buildType", "-DCMAKE_INSTALL_PREFIX=$Install",
    "-DCMAKE_C_COMPILER=$(Join-Path $bin 'cc.exe')", "-DCMAKE_CXX_COMPILER=$(Join-Path $bin 'c++.exe')",
    "-DCMAKE_Fortran_COMPILER=$(Join-Path $bin 'gfortran.exe')",
    "-DMPI_C_COMPILER=$(Join-Path $bin 'mpicc.exe')", "-DMPI_CXX_COMPILER=$(Join-Path $bin 'mpicxx.exe')", "-DMPI_Fortran_COMPILER=$(Join-Path $bin 'mpif90.exe')",
    '-DWITH_MPI=ON', '-DWITH_Hypre=ON', '-DWITH_Mumps=ON', '-DCPACK_BUNDLE_EXTRA_WINDOWS_DLLS=OFF', "-DHypre_INCLUDE_DIR=$include",
    "-DHypre_LIBRARIES=$(Join-Path $lib 'libHYPRE.dll.a')", "-DMumps_INCLUDE_DIR=$include",
    '-DMumps_LIBRARIES=-lmumps-dmo -lmumps-zmo -lmumps-smo -lmumps-cmo', "-DBLAS_openblas_LIBRARY=$openBlas",
    "-DLAPACK_LIB=$openBlas", "-DCMAKE_Fortran_FLAGS=$fortranFlags"
)
$oldPath = $env:PATH
try {
    $env:PATH = "$bin;$oldPath"
    & $cmake @configureArgs
    if ($LASTEXITCODE) { exit $LASTEXITCODE }
    & $cmake --build $Build --target install --parallel
    exit $LASTEXITCODE
} finally { $env:PATH = $oldPath }
