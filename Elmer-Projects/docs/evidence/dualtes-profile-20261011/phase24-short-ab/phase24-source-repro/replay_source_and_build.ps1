param(
  [Parameter(Mandatory=$true)][string]$ElmerSource,
  [Parameter(Mandatory=$true)][string]$BuildDir,
  [Parameter(Mandatory=$true)][string]$InstallPrefix,
  [string]$ToolchainBin = 'G:\Elmer-Migration\msys64\ucrt64\bin',
  [string]$ElmerRuntimeRoot = 'G:\Elmer-Migration\elmer',
  [string]$MpiExec = 'C:\Program Files\Microsoft MPI\Bin\mpiexec.exe',
  [switch]$Build
)
$ErrorActionPreference = 'Stop'
$base = '5a8de867068be0568f09af40fb90fee300dfbede'
$patch = Join-Path $PSScriptRoot '01_phase24_source_from_release_26_1.patch'
$manifest = Get-Content (Join-Path $PSScriptRoot 'source_hash_manifest.json') -Raw | ConvertFrom-Json
$actualPatchHash = (Get-FileHash -Algorithm SHA256 $patch).Hash.ToLowerInvariant()
if ($actualPatchHash -ne $manifest.patch.sha256) { throw "Patch SHA-256 mismatch: $actualPatchHash" }
$head = (git -C $ElmerSource rev-parse HEAD).Trim()
if ($head -ne $base) { throw "Expected clean Elmer $base; got $head" }
$status = git -C $ElmerSource status --porcelain
if ($status) { throw 'ElmerSource must be clean before applying the patch.' }
foreach ($file in $manifest.files) {
  $blob = (git -C $ElmerSource rev-parse ("{0}:{1}" -f $base,$file.path)).Trim()
  if ($blob -ne $file.before_git_blob) { throw "Base hash mismatch: $($file.path)" }
}
git -C $ElmerSource apply --check $patch
if ($LASTEXITCODE -ne 0) { throw 'git apply --check failed; no patch was applied.' }
git -C $ElmerSource apply $patch
if ($LASTEXITCODE -ne 0) { throw 'git apply failed.' }
foreach ($file in $manifest.files) {
  $actual = (git -C $ElmerSource hash-object ("--path={0}" -f $file.path) (Join-Path $ElmerSource $file.path)).Trim()
  if ($actual -ne $file.target_git_blob) { throw "Post-apply hash mismatch: $($file.path)" }
}
Write-Output "SOURCE_REPLAY_PASS files=$($manifest.files.Count)"
if (-not $Build) { return }

$msysRoot = Split-Path $ToolchainBin -Parent
$env:PATH = "$ToolchainBin;$ElmerRuntimeRoot\bin;C:\Program Files\Microsoft MPI\Bin;$env:PATH"
& cmake -S $ElmerSource -B $BuildDir -G Ninja `
  -DCMAKE_BUILD_TYPE=Release `
  -DCMAKE_C_COMPILER="$ToolchainBin\cc.exe" `
  -DCMAKE_CXX_COMPILER="$ToolchainBin\c++.exe" `
  -DCMAKE_Fortran_COMPILER="$ToolchainBin\gfortran.exe" `
  -DCMAKE_MAKE_PROGRAM="$ToolchainBin\ninja.exe" `
  -DCMAKE_INSTALL_PREFIX=$InstallPrefix `
  -DELMER_SOLVER_HOME="$InstallPrefix/share/elmersolver" `
  -DELMER_INSTALL_LIB_DIR=lib/elmersolver `
  -DWITH_MPI=ON `
  -DMPI_C_COMPILER="$ToolchainBin\mpicc.exe" `
  -DMPI_CXX_COMPILER="$ToolchainBin\mpicxx.exe" `
  -DMPI_Fortran_COMPILER="$ToolchainBin\mpif90.exe" `
  -DMPIEXEC_EXECUTABLE=$MpiExec -DMPIEXEC_NUMPROC_FLAG=-n -DMPIEXEC_MAX_NUMPROCS=6 `
  -DWITH_Hypre=ON `
  -DHypre_INCLUDE_DIR="$msysRoot\include" `
  -DHypre_LIBRARIES="$msysRoot\lib\libHYPRE.dll.a" `
  -DWITH_OpenMP=ON -DWITH_LUA=ON -DWITH_UMFPACK=ON -DWITH_Mumps=OFF `
  -DWITH_ELMERGUI=OFF -DWITH_ELMERPOST=OFF -DWITH_ELMERGUILOGGER=OFF `
  -DWITH_ELMERGUITESTER=OFF -DWITH_CONTRIB=OFF `
  -DBLAS_openblas_LIBRARY="$ElmerRuntimeRoot\bin\libopenblas.dll.a"
if ($LASTEXITCODE -ne 0) { throw 'CMake configure failed.' }
& cmake --build $BuildDir --parallel 6
if ($LASTEXITCODE -ne 0) { throw 'Release build failed.' }
Write-Output 'RELEASE_BUILD_PASS'
