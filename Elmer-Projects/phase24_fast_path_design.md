# Phase24 CPU fast path

The Phase24 path is opt-in through `Phase24 Vector Assembly = Logical True` and is limited to a scalar, three-dimensional, fixed-mesh P1 tetrahedral bulk mesh. Unsupported element types, dimensions, DOF counts, periodic flips, anisotropic or MATC/UDF properties, dynamic heat capacity, body-force elements, convection, phase change, and other ineligible terms remain on the original HeatSolve path.

## Persistent initialization state

`HeatSolve.F90` builds contiguous per-bulk-element metadata once per mesh:

- `Phase24ElementNodes(4,nelem)` and `Phase24ElementEq(4,nelem)`
- body/part identifiers
- physical shape gradients and tetra volume
- reusable static stiffness and spatial mass matrices
- compact `StaticElements`, `DynamicElements`, and `GenericElements` lists

`ElementInfo` is called during metadata construction only. The generic `ElementDescription` implementation is unchanged; this avoids affecting periodic-projector and boundary code.

Static elements require constant scalar density, heat capacity, and isotropic conductivity, with no source or temperature-dependent term. Their conductivity stiffness and spatial mass are computed once and reused. The timestep/BDF multiplier and RHS history remain per-iteration quantities. Dynamic and generic elements continue through the existing local formulation.

## Assembly and insertion

Eligible static elements use an OpenMP `parallel do`. Each worker forms private local mass, stiffness, and RHS values. CSR destinations are cached in `CRSMatrix.F90`; the cold cache build is narrowly synchronized, while subsequent matrix updates use cached positions and OpenMP atomics. No whole-element loop is enclosed by a global critical region. A coloring pass is not required for this stage because the atomic implementation is the simpler safe option and removes repeated CRS searches.

The cache records CSR row/column storage identity and dimensions. It is cleared when the matrix structure or mesh changes. The observed one-step run built 266,569 element destinations once, reached 533,138 hits after the second nonlinear iteration, and reported only the two expected initial invalidations.

## Fallback and validation

The fallback path is explicitly guarded against Fortran non-short-circuit evaluation: fast-path arrays are never indexed until activation and bounds have been checked. The production hybrid prism/tetra case therefore disables the fast path cleanly and completes assembly through the generic route.

The reported conformal CPU HYPRE smoke uses `case_dual_pulse_center` on the tetra-only `mesh_dual_base` mesh. One timestep converged with FlexGMRES + BoomerAMG. Solver wall time improved from 20.75 s to 11.18 s; assembly improved from 13.51 s to 2.23 s. Matrix sparsity was identical, normalized matrix-value error was `2.49e-8`, RHS error was `1.17e-16`, and the maximum temperature-field difference was `1.55e-4` in the recorded one-step outputs.

## Build note

The existing configured native build used MSYS2 UCRT64 GNU Fortran/C/C++ wrappers with Microsoft MPI import libraries and `mpiexec` from Microsoft MPI. The CMake cache already contained the correct MPI paths (`C:/msys64/ucrt64/.../mpif90.exe`, `mpicc.exe`, `mpicxx.exe`, `libmsmpi.dll.a`, and `C:/Program Files/Microsoft MPI/Bin/mpiexec.exe`). The apparent compiler failure was the MSYS2 runtime missing from the PowerShell `PATH`: `f951.exe` exited with Windows status `0xC0000139`. Prepending `C:\msys64\ucrt64\bin;C:\msys64\usr\bin` restored the configured build without reconfiguration or MPI changes. Mixing the WSL OpenMPI wrappers with the native cache was avoided.

GPU assembly was not changed.
