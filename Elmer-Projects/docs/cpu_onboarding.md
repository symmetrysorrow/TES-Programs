# CPU-only setup and generated-mesh pilot (experimental)

This is the first external-user onboarding path. It does **not** claim that
a freshly generated mesh reproduces the validated fixed `h8` single-pixel
waveform. Do not publish physical results from the generated-mesh pilot until
mesh convergence and COMSOL/reference parity have been checked.

## Supported setup (first milestone)

- Windows 11 + WSL2, Ubuntu 24.04 x86_64 (WSL distro name `ElmerGPU` by default)
- Windows Python 3.11+ (the Python standard-library `venv` module is required)
- Network access, sufficient disk space and permission to run `apt-get` in WSL
- CPU-only solver; **no** AMD/NVIDIA GPU, ROCm, HIP or CUDA required

The setup command currently expects Ubuntu 24.04 already installed in WSL.
It does not install Windows, Python, WSL itself or register/rename the distro.
If the distro has another name, set `TES_WSL_DISTRO` in the terminal before
running setup or calculation.

From the `Elmer-Projects` directory on Windows:

```powershell
python main.py setup
.\.venv\Scripts\python.exe main.py list
.\.venv\Scripts\python.exe main.py run single_pixel_generated --dry-run
.\.venv\Scripts\python.exe main.py mesh single_pixel_generated
.\.venv\Scripts\python.exe main.py run single_pixel_generated
```

`setup` creates a host Python virtual environment, downloads the declared
Python requirements, then invokes `src/support/setup_cpu_ubuntu.sh` as
root inside WSL. The Ubuntu installer retrieves the pinned Elmer/HYPRE
source commits, verifies and applies the repository patches, installs Ubuntu
build dependencies, compiles a **CPU-only** solver, and builds the TES
Fortran user function. Build outputs are isolated under
`/opt/tes-elmer-cpu`; source and build trees are separate. Ubuntu and Python
package versions actually installed are saved alongside the build.

The source lock and changes are in
[`elmergpu-sources.lock`](elmergpu-sources.lock),
[`elmer_phase24_custom.patch`](elmer_phase24_custom.patch) and
[`hypre_hip_rocm714_compat.patch`](hypre_hip_rocm714_compat.patch).
The HYPRE patch is reused even for the CPU build, but the HIP and CUDA
compiler options are explicitly disabled.

**Caveats**: This installation process has not yet been verified on a fresh
third-party PC. Package versions are not fully locked; setup may download
several gigabytes and requires elevated privileges. Neither `setup`
nor the Python tests establish numerical parity. Check build logs, produced
SIF files, runtime provenance, and published reference data before relying on
results.

## Mesh settings

The existing `single_pixel_alpha240` project retains the validated,
prebuilt `mesh_hybrid_fullconf_h8` reference. That mesh is not distributed
with this repository and must not be silently replaced.

The new opt-in `single_pixel_generated` project points to
`models/single_pixel_generated.toml` and has **no** `mesh.source` field.
Its `[mesh]` table specifies global/local sizes and layer counts. Running
`mesh` or `run` generates Gmsh geometry and converts it with ElmerGrid.
On a freshly bootstrapped Windows host, if no native ElmerGrid exists, the
converter is taken from the bootstrapped WSL CPU installation.

A previously generated mesh is only reused when its recorded TOML
mesh hash matches. Missing or inconsistent provenance causes a fail-closed
error; the user must move the old cache aside to regenerate it. Generated
meshes and results remain out of Git.

## Ubuntu directly

On Ubuntu 24.04 x86_64, the **solver build only** can be started via:

```bash
sudo bash src/support/setup_cpu_ubuntu.sh
```

The user-facing `main.py run` path still targets WSL from Windows.
Native Linux execution remains a separate integration task.

## Remaining release gates

1. Verify clean-machine CPU build and run end-to-end, including `ElmerGrid`.
2. Verify `single_pixel_generated` connectivity, mesh-convergence and
   waveform parity against fixed h8/COMSOL data.
3. Extend the new bootstrap to auto-provision WSL or to execute natively
   on Linux; publish a reproducible package lock and acceptance tests.
