[CmdletBinding()]
param(
    [ValidateSet("cuda", "hip")]
    [string]$Backend = "cuda",
    [string]$GpuArchitecture = "86",
    [string]$HypreTag = "v3.0.0",
    [switch]$BuildExamples
)

$ErrorActionPreference = "Stop"

$repo = (Resolve-Path (Join-Path $PSScriptRoot "..\..")).Path
$tools = Join-Path (Split-Path $repo -Parent) "tools"
if ($repo -notmatch '^[A-Za-z]:\\') { throw "Expected a repository on a Windows drive" }
$drive = $repo.Substring(0, 1).ToLowerInvariant()
$repoWsl = "/mnt/$drive" + (($repo.Substring(2) -replace '\\', '/'))
$toolsWsl = "/mnt/$drive" + (($tools.Substring(2) -replace '\\', '/'))

$hypreSrc = "$toolsWsl/hypre-$Backend-src"
$hypreBuild = "/home/symme/hypre-$Backend-build"
$hypreInstall = "$toolsWsl/hypre-$Backend-install"
$elmerBuild = "/home/symme/elmer-hypre-$Backend-build"
$elmerInstall = "$toolsWsl/elmer-hypre-$Backend-wsl"
$amgxInstall = "$toolsWsl/amgx-gpu-install-mpi"
$gpuOptions = if ($Backend -eq "cuda") {
    "-DHYPRE_ENABLE_CUDA=ON -DCMAKE_CUDA_ARCHITECTURES=$GpuArchitecture -DCMAKE_CUDA_HOST_COMPILER=/usr/bin/g++-13 -DHYPRE_ENABLE_CUSPARSE=ON"
} else {
    "-DHYPRE_ENABLE_HIP=ON -DCMAKE_HIP_ARCHITECTURES=$GpuArchitecture"
}
$examples = if ($BuildExamples) { "ON" } else { "OFF" }

# Do not touch the existing CPU or AMGX prefixes.  This creates a separately
# reproducible backend-specific HYPRE + Elmer installation.
$bash = @"
set -euo pipefail
# CMake otherwise searches the inherited, very long Windows PATH once per CUDA
# component; retain only the WSL toolchains required by either backend.
export PATH=/opt/rocm-wsl/bin:/usr/lib/nvidia-cuda-toolkit/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/usr/lib/wsl/lib
if [ '$Backend' = cuda ]; then
  # CUDA 12.4 in this WSL image supports GCC <= 13, while /usr/bin/c++ is GCC 15.
  export CUDAHOSTCXX=/usr/bin/g++-13
fi
if [ ! -d '$hypreSrc/.git' ]; then
  git clone --depth 1 --branch '$HypreTag' https://github.com/hypre-space/hypre.git '$hypreSrc'
fi
cmake -S '$hypreSrc/src' -B '$hypreBuild' -G Ninja \
  -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX='$hypreInstall' \
  -DBUILD_SHARED_LIBS=ON \
  -DHYPRE_ENABLE_MPI=ON -DHYPRE_BUILD_EXAMPLES=$examples \
  -DHYPRE_ENABLE_UNIFIED_MEMORY=OFF -DHYPRE_ENABLE_GPU_AWARE_MPI=OFF \
  -DHYPRE_ENABLE_UMPIRE=OFF $gpuOptions
cmake --build '$hypreBuild' --parallel
cmake --install '$hypreBuild'
cmake -S '$toolsWsl/elmer-hypre/src' -B '$elmerBuild' -G Ninja \
  -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX='$elmerInstall' \
  -DWITH_MPI=ON -DWITH_Mumps=ON -DWITH_Hypre=ON -DWITH_AMGX=ON \
  -DMUMPSROOT=/usr -DPARMETISROOT=/usr \
  -DMetis_INCLUDE_DIR=/usr/include -DMetis_LIBRARIES=/usr/lib/x86_64-linux-gnu/libmetis.so \
  -DParMetis_INCLUDE_DIR=/usr/include/parmetis -DParMetis_LIBRARIES=/usr/lib/x86_64-linux-gnu/libparmetis.so \
  -DHypre_INCLUDE_DIR='$hypreInstall/include' -DHypre_LIBRARIES='$hypreInstall/lib/libHYPRE.so' \
  -DAMGXINCLUDE='$amgxInstall/include' -DAMGXLIB='$amgxInstall/lib'
cmake --build '$elmerBuild' --parallel
cmake --install '$elmerBuild'
"@

Write-Host "Building the portable HYPRE $Backend backend in WSL."
& wsl.exe -d Ubuntu -- bash -lc $bash
if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }
Write-Host "Installed: $elmerInstall"
