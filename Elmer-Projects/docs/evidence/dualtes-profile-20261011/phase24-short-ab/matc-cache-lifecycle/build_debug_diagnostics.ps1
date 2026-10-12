param([Parameter(Mandatory=$true)][string]$ElmerSource,[Parameter(Mandatory=$true)][string]$BuildDir,[Parameter(Mandatory=$true)][string]$InstallPrefix,[string]$ToolchainBin='G:\Elmer-Migration\msys64\ucrt64\bin',[string]$ElmerRuntimeRoot='G:\Elmer-Migration\elmer',[string]$MpiExec='C:\Program Files\Microsoft MPI\Bin\mpiexec.exe')
$ErrorActionPreference='Stop'
$PSNativeCommandUseErrorActionPreference=$false
$msysRoot=Split-Path $ToolchainBin -Parent
$env:PATH="$ToolchainBin;$ElmerRuntimeRoot\bin;C:\Program Files\Microsoft MPI\Bin;$env:PATH"
$env:CMAKE_Fortran_FLAGS_DEBUG='-O0 -g -fcheck=all -finit-real=snan -finit-integer=-999 -fbacktrace -ffpe-trap=invalid,zero,overflow'
$env:CMAKE_C_FLAGS_DEBUG='-O0 -g -fno-omit-frame-pointer'
$env:CMAKE_CXX_FLAGS_DEBUG='-O0 -g -fno-omit-frame-pointer'
$ErrorActionPreference='Continue'
& cmake -S $ElmerSource -B $BuildDir -G Ninja `
  -DCMAKE_BUILD_TYPE=Debug -DCMAKE_C_COMPILER="$ToolchainBin\cc.exe" `
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
if($LASTEXITCODE){throw 'Debug configure failed'}
$ErrorActionPreference='Stop'
Write-Output "DEBUG_FLAGS=$env:CMAKE_Fortran_FLAGS_DEBUG"
$ErrorActionPreference='Continue'
& cmake --build $BuildDir --target Solver_TGT HeatSolve --parallel 6
if($LASTEXITCODE){throw 'Debug solver/HeatSolve build failed'}
$ErrorActionPreference='Stop'
Write-Output 'DEBUG_BOUNDS_CHECK_BUILD_PASS'

