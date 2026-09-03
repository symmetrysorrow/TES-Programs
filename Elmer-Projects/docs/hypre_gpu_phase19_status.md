# HYPRE CUDA Phase19 status

## Implemented

- Added generated FlexGMRES + BoomerAMG CPU/CUDA/HIP case definitions.
- Added one-step matrix/RHS diagnostics and an exact record comparator.
- Added a WSL CUDA/HIP build launcher and runtime launcher.
- Added `HYPRE GPU` input registration to the local Elmer build and selected
  HYPRE device execution before IJ construction.  The runtime sets
  `CUDA_VISIBLE_DEVICES=0` (or `HIP_VISIBLE_DEVICES=0`) for the one-rank test.
- FlexGMRES now aborts when HYPRE reaches its iteration limit without meeting
  the requested residual, instead of returning a silently inaccurate iterate.
- Added CPU and CUDA one-step no-mortar controls.  These controls use the same
  restart, mesh, source, and solver tolerances but disable the explicit mortar
  rows, so they test the ordinary SPD HYPRE path independently.

The Elmer source tree under `tools/` is ignored by the project repository; the
local build used the corresponding edits in `tools/elmer-hypre/src/fem/src/`
(`SolveHypre.c`, `SParIterSolver.F90`, and `SOLVER.KEYWORDS`).  Those edits
must be carried into the Elmer source tree or an upstream patch when packaging
the change.

## Reproducible checks

The one-step CPU and CUDA runs both loaded the same restart and assembled
1,447,196 matrix records and 87,534 RHS records.  The A and b dumps matched
exactly (record count, mismatch count, maximum absolute difference, and
SHA-256).

The GPU path initialized the RTX 3060 Ti successfully.  On the mortar system,
however, FlexGMRES + BoomerAMG reached 2,000 iterations with relative residual
`1.81166e-4` (requested `1e-11`), so the fail-fast guard correctly terminated
the case.  The prior run without the guard produced a different TES result
from MUMPS despite identical A/b; that result must not be accepted.

This confirms that the remaining issue is the preconditioner/formulation for
the mortar-constrained system, not CUDA migration or matrix assembly.  An
explicit saddle-point formulation requires an MGR-style block treatment (or a
validated exact elimination/reduced operator) before full transient timing.

## Correction policy and re-check

1. Treat the assembled A,b pair as the invariant.  CPU/GPU dumps must match
   exactly before comparing solver results.
2. Treat a constrained solve as successful only when the reported relative
   residual is at or below `1e-11`; iteration-limit exits are failures.
3. Do not enable `Eliminate Linear Constraints` as a proxy for elimination:
   this case reports `Number of eliminated rows: 0`, so the matrix remains an
   explicit saddle system.
4. Keep the mortar case disabled for BoomerAMG production until C/F row
   markers are passed into an MGR block preconditioner (or an exact reduced
   operator is validated against MUMPS).

The new no-mortar controls validate this boundary: both CPU and CUDA runs
converged in 286 FlexGMRES iterations to relative residual `9.630026e-12`,
with identical result norm `0.15677291527764314`.  Their A dumps contain
1,340,068 records and their b dumps 84,636 records; CPU and GPU hashes are
identical (`e8ae9210…` and `cbb51c2c…`, respectively).  The CUDA log also
confirms execution on an NVIDIA GeForce RTX 3060 Ti.

The final constrained CPU re-check (with no elimination keyword) exits 1 as
intended: it enters `SolveWithLinearRestriction` and FlexGMRES stops at 2,000
iterations with relative residual `2.16316e-7`, still above `1e-11`.  The
earlier GPU run with the ineffective elimination request stopped at
`1.81166e-4`; both are failures.  No constrained result is therefore accepted
as a solution or as a performance measurement.

## MGR integration (option 1)

The local Elmer/HYPRE bridge now exposes HYPRE MGR as preconditioner index 12
(reported by Elmer as method index 912 when combined with FlexGMRES).  For the
explicit mortar system, `SolveWithLinearRestriction` passes the contiguous
constraint-row range to the C layer, which builds the MGR C/F marker array with
`HYPRE_MGRSetCpointsByPointMarkerArray` and keeps the primal rows as F-points.

The one-rank CPU and CUDA MGR smoke cases both initialize successfully.  They
report C rows `84637..87534` (2,898 rows), and the CUDA case initializes the
RTX 3060 Ti.  Adding a dedicated BoomerAMG F solver and changing restriction
to approximate-inverse (`R=3`, `P=2`) improves both CPU and CUDA to the same
relative residual `2.12417e-10` after 2,000 iterations (requested `1e-11`).
The fail-fast guard returns exit code 1, so no MGR solution is accepted.

Thus option 1 is implemented and the library path is exercised on CPU/GPU,
but numerical verification is not yet passed.  The next required work is
calibration of the MGR F-relaxation/coarse (Schur) solver for this mortar
projector, followed by an MPI-rank validation; loosening the tolerance or
silently accepting the final iterate would be incorrect.

## Block diagnosis and tuning results

The saved one-step matrix was split at row 84,637 (`n_F=84,636`,
`n_C=2,898`).  The measured blocks are:

| quantity | value |
|---|---:|
| `nnz(K), nnz(B), nnz(Bt), nnz(D)` | `1,340,068`, `53,564`, `53,564`, `0` |
| `||D||` | `0` |
| `||B-Bt^T||/||B||` | `0` |
| `||K-K^T||/||K||` | `1.0400e-5` |
| `diag(K)` min/max | `1.61847e-10 / 6.166998e-2` |
| K row norm min/max | `1.64079e-10 / 9.32446e-2` |
| B row norm min/max | `5.61702e-11 / 2.14576e-10` |

The `K` block is byte-for-byte identical to the no-mortar matrix dump.  This
rules out the F block and assembly as the primary cause; the difficult part is
the very small, zero-diagonal constraint Schur block
`S = -B K^{-1} B^T`.

An independent one-rank block-Schur reference (SuperLU factorization of `K`,
dense `S`) gives `||Ax-b||/||b|| = 3.36e-11` and absolute constraint residual
`6.70e-25`; the reconstructed primal norm is `45.75255699` versus the saved
MUMPS primal norm `45.75255316` (relative difference `2.21e-6`, attributable to
the ill-conditioned Schur reduction).  The calculation is recorded in
`results/block_schur_reference.json`.

MGR tuning was then limited to targeted comparisons.  Registering a dedicated
GPU-capable BoomerAMG F solver reduced the baseline residual from `4.25e-6` to
`5.68e-10`.  Approximate-inverse restriction (`R=3`) improved it further to
`2.12e-10`, while stronger F cycles, approximate-inverse interpolation
(`P=4`), coarse-solver replacement, and symmetric constraint scaling all
degraded the result.  The best tested one-level configuration is therefore
dedicated BoomerAMG F (`tol=1e-4`, max 2 AMG iterations), `R=3`, `P=2`, with
`NonCpointsToFpoints=1`; it still misses the required `1e-11` by about 21x.

HYPRE 3.1.0 and MPI 2/4-rank validation have not yet been promoted to the
baseline.  The current conclusion is that MGR is a valid library/GPU path but
is not yet a production solver for this mortar system.  The next production
candidate is an explicit Schur preconditioner: BoomerAMG for `K` (GPU) and a
small, separately solved constraint Schur problem.
