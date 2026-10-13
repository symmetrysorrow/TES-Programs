# Rebuilding the Elmer solver in WSL (ElmerGPU distro)

The HIP build uses upstream source revisions pinned in
[`elmergpu-sources.lock`](elmergpu-sources.lock), with TES/Phase24 changes and
the HYPRE HIP compatibility fix carried as reviewable patches in this repo.
Build outputs stay in the WSL distro and are not committed.

## Host requirements

- Windows 11 with WSL2 and an AMD GPU supported by AMD's ROCm WSL stack.
- An Ubuntu 24.04 WSL distribution. The existing runner expects the distro
  name `ElmerGPU` and installs under `/opt/elmer-gpu`.
- Ubuntu's `multiverse` apt component enabled for the pinned `libparmetis-dev`
  package.
- ROCm/HIP 7.2.4 with the matching `librocdxg` component. The build currently
  targets `gfx1201`; another GPU architecture needs a reviewed lock/build
  change.

Install Ubuntu build packages and validate the ROCm version from the repo root:

```powershell
wsl -d ElmerGPU -u root -- bash /mnt/d/github/TES-Programs/Elmer-Projects/src/support/install_elmergpu_build_deps.sh
wsl -d ElmerGPU -u root -- bash /mnt/d/github/TES-Programs/Elmer-Projects/src/support/build_elmer_hypre_hip_elmergpu.sh
```

Adjust the mounted repo path in these commands if the checkout is elsewhere.
The source preparation helper clones both upstream projects at their locked
base commits and applies the tracked patches under `/root/elmergpu-src`. If you
already have a source directory there with local changes, set
`ELMERGPU_SOURCE_ROOT` to a new empty directory when invoking the build script.

Ubuntu package versions are pinned in
[`elmergpu-ubuntu-24.04-packages.txt`](elmergpu-ubuntu-24.04-packages.txt).
The AMD GPU package versions are recorded in
[`elmergpu-rocm-7.2.4-packages.txt`](elmergpu-rocm-7.2.4-packages.txt); install
the AMD ROCm WSL apt repository first so those pinned packages are available.
The lock list includes `rocdxg-roct`, the WSL GPU bridge runtime. Package
archives remain hosted by Ubuntu and AMD; the repo stores the exact package
names and versions, while the upstream repositories provide the binaries.

## Source changes

- [`elmer_phase24_custom.patch`](elmer_phase24_custom.patch) records all Elmer
  changes relative to upstream base `5a8de867068be0568f09af40fb90fee300dfbede`.
- [`hypre_hip_rocm714_compat.patch`](hypre_hip_rocm714_compat.patch) records
  the HIP/Thrust compatibility changes relative to HYPRE `9dc9e18...`.
