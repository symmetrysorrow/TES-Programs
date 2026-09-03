# Phase 19 HYPRE GPU / mortar status

This is the implementation handoff for the Phase 19 saddle-point experiment.
The acceptance gate remains the original full-system residual tolerance of
`1e-11`; an iteration-limit result is not accepted.

## Verified algebra

The saved explicit mortar matrix has `nF = 84636`, `nC = 2898`,
`nnz(K) = 1,340,068`, `nnz(B) = nnz(B^T) = 53,564`, and `D = 0`.
The block-driver K dump matches the no-mortar K dump byte-for-byte for the
tested case.  The independent SuperLU Schur oracle gives
`||A x - b||_2 = 1.14e-16` and `||B u||_2 = 1.16e-24`.

The saved MUMPS reference is primal-only (`84636` entries), so it is used for
primal agreement, not for a fictitious full saddle residual.  The exact Schur
primal agrees with it to `2.21e-6` relative.

## Candidate results

| Candidate | Result |
|---|---|
| HYPRE 3.0 FlexGMRES + BoomerAMG, no mortar, CPU | converged in 286 iterations, relative residual `9.63e-12` |
| HYPRE 3.0 FlexGMRES + BoomerAMG, no mortar, CUDA-enabled library with GPU disabled/enabled | both matched the 286-iteration CPU result in the existing A/B run |
| HYPRE 3.0 MGR, explicit mortar, CPU/GPU | both hit 2000 iterations at `2.12417e-10`; rejected |
| HYPRE 3.1 MGR, explicit mortar, CPU | NaN from iteration 2; HYPRE reported input INF/NaN and exited 1; rejected |
| HYPRE 3.1 BoomerAMG, no mortar, CPU | converged in 286 iterations, relative residual `9.63e-12`; control passed |
| HYPRE 3.1 BoomerAMG, no mortar, GPU | CUDA runtime reported no CUDA-capable device in this process; not a solver result |
| Native block, diagonal Schur, one-cycle AMG K + direct Schur | initially decreases, then diverges; rejected |
| Native block, Gauss–Seidel variant | diverges earlier; rejected |
| Native block, diagonal Schur, ten-cycle AMG K | still diverges; rejected |

The exact-versus-approximate Schur diagnostic is in
`results/schur_strategy_comparison.json`, with the ILU candidate in
`results/schur_strategy_ilu.json`.  The exact setup is about 173 s and
passes the algebraic gate; diagonal setup is about 1 s but leaves a full
relative residual of about `1.99`; ILU setup is about 82 s and leaves a full
relative residual of about `3.06`.  These are independent CPU reference
experiments; the ILU case is a portable proxy for an approximate AMG action,
not a claim that HYPRE used ILU.

## Implementation

The block candidate is opt-in.  It sets one BoomerAMG cycle for the large
primal block and routes the small Schur approximation to UMFPACK through the
nested solver's own parameter list.  The ignored Elmer source tree is carried
as [hypre_gpu_phase19_blocksolve.patch](hypre_gpu_phase19_blocksolve.patch).

The HYPRE build helper now keeps non-default tags side-by-side.  HYPRE
`v3.1.0` CUDA and Elmer builds complete after explicitly pointing CMake at the
WSL `/usr/include` CUDA headers.  The 3.1 no-mortar CPU control is unchanged,
but 3.1 MGR becomes NaN immediately, so upgrading HYPRE alone is not a fix.
The current WSL CUDA process also reports no CUDA-capable device even though
`nvidia-smi` exposes an RTX 3060 Ti; GPU promotion therefore remains blocked
by the runtime/device mismatch.

The HIP helper now prefers ROCm's Thrust headers over the CUDA Thrust headers
present in `/usr/include`.  The first HIP build exposed a mixed-header failure;
after fixing the include order, HYPRE 3.1 still stops in
`csr_spgemm_device_symbl.c` at `HYPRE_THRUST_IDENTITY(char)` under ROCm 7.14.
Thus HIP is toolchain-blocked, not solver-validated.  ROCm sees an AMD Radeon
RX 9070 XT in this WSL image, so this is an actionable HYPRE/ROCm compatibility
fix for a subsequent pass.

No short transient is accepted yet: the available seven-step run was
interrupted during the rejected MGR investigation.  MPI is not yet a result
for this case because the mesh has no `partitioning.2` directory.  These are
explicitly open validation gates, not successful production evidence.

## Recommendation at this point

Do not promote the 3.0/3.1 MGR result or the native block approximation to
production.  The exact Schur path is a correctness oracle, not a production
memory strategy.  Promotion requires a matrix-free Schur/AMG variant to pass
the full residual and constraint gates, followed by the short-transient, MPI,
and real GPU checks.
