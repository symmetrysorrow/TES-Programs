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
$phaseDir = Join-Path $PSScriptRoot '..\phase24-source-repro'
$patch01 = Join-Path $phaseDir '01_phase24_source_from_release_26_1.patch'
$manifest01 = Get-Content (Join-Path $phaseDir 'source_hash_manifest.json') -Raw | ConvertFrom-Json
$patch02 = Join-Path $PSScriptRoot '02_matc_ast_purity_validator.patch'
$manifest02 = Get-Content (Join-Path $PSScriptRoot '02_source_hash_manifest.json') -Raw | ConvertFrom-Json
if ((Get-FileHash $patch01 -Algorithm SHA256).Hash.ToLowerInvariant() -ne $manifest01.patch.sha256) { throw '01 patch hash mismatch.' }
if ((Get-FileHash $patch02 -Algorithm SHA256).Hash.ToLowerInvariant() -ne $manifest02.patch.sha256) { throw '02 patch hash mismatch.' }
if ((git -C $ElmerSource rev-parse HEAD).Trim() -ne $base) { throw "Expected upstream $base." }
if (git -C $ElmerSource status --porcelain) { throw 'ElmerSource must be clean.' }
foreach ($f in $manifest01.files) {
  if ((git -C $ElmerSource rev-parse ("{0}:{1}" -f $base,$f.path)).Trim() -ne $f.before_git_blob) { throw "01 preimage mismatch: $($f.path)" }
}
git -C $ElmerSource apply --check --whitespace=nowarn $patch01
if ($LASTEXITCODE) { throw '01 git apply --check failed.' }
git -C $ElmerSource apply --whitespace=nowarn $patch01
if ($LASTEXITCODE) { throw '01 apply failed.' }
foreach ($f in $manifest01.files) {
  if ((git -C $ElmerSource hash-object ("--path={0}" -f $f.path) (Join-Path $ElmerSource $f.path)).Trim() -ne $f.target_git_blob) { throw "01 postimage mismatch: $($f.path)" }
}
Write-Output "PHASE24_PATCH_01_PASS files=$($manifest01.files.Count)"
foreach ($f in $manifest02.files) {
  if ((git -C $ElmerSource hash-object ("--path={0}" -f $f.path) (Join-Path $ElmerSource $f.path)).Trim() -ne $f.before_git_blob) { throw "02 preimage mismatch: $($f.path)" }
}
git -C $ElmerSource apply --check --whitespace=nowarn $patch02
if ($LASTEXITCODE) { throw '02 git apply --check failed.' }
git -C $ElmerSource apply --whitespace=nowarn $patch02
if ($LASTEXITCODE) { throw '02 apply failed.' }
foreach ($f in $manifest02.files) {
  if ((git -C $ElmerSource hash-object ("--path={0}" -f $f.path) (Join-Path $ElmerSource $f.path)).Trim() -ne $f.after_git_blob) { throw "02 postimage mismatch: $($f.path)" }
}
Write-Output "MATC_AST_PATCH_02_PASS files=$($manifest02.files.Count)"
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
if ($LASTEXITCODE) { throw 'CMake configure failed.' }
& cmake --build $BuildDir --clean-first --parallel 6
if ($LASTEXITCODE) { throw 'Release build failed.' }
Write-Output 'AST_RELEASE_BUILD_PASS'
