param(
  [Parameter(Mandatory=$true)][string]$ElmerSource,
  [Parameter(Mandatory=$true)][string]$BuildDir,
  [Parameter(Mandatory=$true)][string]$InstallPrefix,
  [string]$ToolchainBin='G:\Elmer-Migration\msys64\ucrt64\bin',
  [string]$ElmerRuntimeRoot='G:\Elmer-Migration\elmer',
  [string]$MpiExec='C:\Program Files\Microsoft MPI\Bin\mpiexec.exe',
  [switch]$Build
)
$ErrorActionPreference='Stop'
$PSNativeCommandUseErrorActionPreference=$false
$base='5a8de867068be0568f09af40fb90fee300dfbede'
$phase=Join-Path $PSScriptRoot '..\phase24-source-repro'
$ast=Join-Path $PSScriptRoot '..\matc-ast-release'
$p1=Join-Path $phase '01_phase24_source_from_release_26_1.patch'
$m1=Get-Content (Join-Path $phase 'source_hash_manifest.json') -Raw | ConvertFrom-Json
$p2=Join-Path $ast '02_matc_ast_purity_validator.patch'
$m2=Get-Content (Join-Path $ast '02_source_hash_manifest.json') -Raw | ConvertFrom-Json
$p3=Join-Path $PSScriptRoot '03_matc_cache_lifecycle.patch'
$m3=Get-Content (Join-Path $PSScriptRoot '03_source_hash_manifest.json') -Raw | ConvertFrom-Json
foreach($x in @(@($p1,$m1.patch.sha256),@($p2,$m2.patch.sha256),@($p3,$m3.patch.sha256))) {
  if ((Get-FileHash $x[0] -Algorithm SHA256).Hash.ToLowerInvariant() -ne $x[1]) { throw "Patch SHA mismatch: $($x[0])" }
}
if ((git -C $ElmerSource rev-parse HEAD).Trim() -ne $base) { throw 'Elmer source is not the pinned official revision.' }
if (git -C $ElmerSource status --porcelain) { throw 'Elmer source must be clean before replay.' }
foreach($f in $m1.files) { if((git -C $ElmerSource rev-parse ("{0}:{1}" -f $base,$f.path)).Trim() -ne $f.before_git_blob) { throw "01 base blob mismatch: $($f.path)" } }
git -C $ElmerSource apply --check --whitespace=nowarn $p1; if($LASTEXITCODE) { throw '01 apply check failed' }
git -C $ElmerSource apply --whitespace=nowarn $p1; if($LASTEXITCODE) { throw '01 apply failed' }
foreach($f in $m1.files) { if((git -C $ElmerSource hash-object ("--path={0}" -f $f.path) (Join-Path $ElmerSource $f.path)).Trim() -ne $f.target_git_blob) { throw "01 postimage mismatch: $($f.path)" } }
foreach($f in $m2.files) { if((git -C $ElmerSource hash-object ("--path={0}" -f $f.path) (Join-Path $ElmerSource $f.path)).Trim() -ne $f.before_git_blob) { throw "02 preimage mismatch: $($f.path)" } }
git -C $ElmerSource apply --check --whitespace=nowarn $p2; if($LASTEXITCODE) { throw '02 apply check failed' }
git -C $ElmerSource apply --whitespace=nowarn $p2; if($LASTEXITCODE) { throw '02 apply failed' }
foreach($f in $m2.files) { if((git -C $ElmerSource hash-object ("--path={0}" -f $f.path) (Join-Path $ElmerSource $f.path)).Trim() -ne $f.after_git_blob) { throw "02 postimage mismatch: $($f.path)" } }
foreach($f in $m3.files) {
  $target=Join-Path $ElmerSource $f.path
  if((git -C $ElmerSource hash-object ("--path={0}" -f $f.path) $target).Trim() -ne $f.before_git_blob) { throw "03 preimage blob mismatch: $($f.path)" }
  if((Get-FileHash $target -Algorithm SHA256).Hash.ToLowerInvariant() -ne $f.before_raw_sha256) { throw "03 preimage SHA mismatch: $($f.path)" }
}
git -C $ElmerSource apply --check --whitespace=nowarn $p3; if($LASTEXITCODE) { throw '03 apply check failed' }
git -C $ElmerSource apply --whitespace=nowarn $p3; if($LASTEXITCODE) { throw '03 apply failed' }
foreach($f in $m3.files) {
  $target=Join-Path $ElmerSource $f.path
  if((git -C $ElmerSource hash-object ("--path={0}" -f $f.path) $target).Trim() -ne $f.after_git_blob) { throw "03 postimage blob mismatch: $($f.path)" }
  if((Get-FileHash $target -Algorithm SHA256).Hash.ToLowerInvariant() -ne $f.after_raw_sha256) { throw "03 postimage SHA mismatch: $($f.path)" }
}
Write-Output 'PATCH_SEQUENCE_01_02_03_PASS'
if(-not $Build) { return }
$msysRoot=Split-Path $ToolchainBin -Parent
$env:PATH="$ToolchainBin;$ElmerRuntimeRoot\bin;C:\Program Files\Microsoft MPI\Bin;$env:PATH"
$ErrorActionPreference='Continue'
& cmake -S $ElmerSource -B $BuildDir -G Ninja `
  -DCMAKE_BUILD_TYPE=Release -DCMAKE_C_COMPILER="$ToolchainBin\cc.exe" `
  -DCMAKE_CXX_COMPILER="$ToolchainBin\c++.exe" -DCMAKE_Fortran_COMPILER="$ToolchainBin\gfortran.exe" `
  -DCMAKE_MAKE_PROGRAM="$ToolchainBin\ninja.exe" -DCMAKE_INSTALL_PREFIX=$InstallPrefix `
  -DELMER_SOLVER_HOME="$InstallPrefix/share/elmersolver" -DELMER_INSTALL_LIB_DIR=lib/elmersolver `
  -DWITH_MPI=ON -DMPI_C_COMPILER="$ToolchainBin\mpicc.exe" -DMPI_CXX_COMPILER="$ToolchainBin\mpicxx.exe" `
  -DMPI_Fortran_COMPILER="$ToolchainBin\mpif90.exe" -DMPIEXEC_EXECUTABLE=$MpiExec `
  -DMPIEXEC_NUMPROC_FLAG=-n -DMPIEXEC_MAX_NUMPROCS=6 -DWITH_Hypre=ON `
  -DHypre_INCLUDE_DIR="$msysRoot\include" -DHypre_LIBRARIES="$msysRoot\lib\libHYPRE.dll.a" `
  -DWITH_OpenMP=ON -DWITH_LUA=ON -DWITH_UMFPACK=ON -DWITH_Mumps=OFF `
  -DWITH_ELMERGUI=OFF -DWITH_ELMERPOST=OFF -DWITH_ELMERGUILOGGER=OFF `
  -DWITH_ELMERGUITESTER=OFF -DWITH_CONTRIB=OFF -DBLAS_openblas_LIBRARY="$ElmerRuntimeRoot\bin\libopenblas.dll.a"
if($LASTEXITCODE) { throw 'CMake configure failed' }
$ErrorActionPreference='Stop'
$ErrorActionPreference='Continue'
& cmake --build $BuildDir --clean-first --parallel 6
if($LASTEXITCODE) { throw 'Release build failed' }
$ErrorActionPreference='Stop'
Write-Output 'PATCH_03_RELEASE_BUILD_PASS'


