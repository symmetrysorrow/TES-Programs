# Phase20 bounded block-Schur measurement status

Updated: 2026-09-05

Phase20 now separates configuration integrity, numerical diagnostics, physical
acceptance, and performance readiness. The four bounded peer cases are lower
and full factorization on CPU and GPU, with the same mesh, restart, one-step
workload, outer limit 15, and Schur settings `tol=1e-4`, `max=30`,
`restart=30`.

| case | backend | factorization | generated timestep | matrix dump |
|---|---|---|---:|---:|
| `case_p20_hypre_block_lower_cpu_probe` | CPU | lower | exactly 1 | off |
| `case_p20_hypre_block_lower_gpu_probe` | GPU | lower | exactly 1 | off |
| `case_p20_hypre_block_full_cpu_probe` | CPU | full | exactly 1 | off |
| `case_p20_hypre_block_full_gpu_probe` | GPU | full | exactly 1 | off |

For Phase20, `matrix_dump_prefix` is retained only as a name. Matrix output
requires the explicit `matrix_dump=true` opt-in; probe cases and the default
parameter sweep set it to false. Legacy Phase19 dump cases retain their
historical prefix behavior for compatibility.

## Measurement contract

The native implementation is in the reviewable patch
`docs/hypre_gpu_phase20_probe.patch`. Its probe is opt-in and resets on a
disabled lifecycle or a changed prefix. It uses monotonic `SYSTEM_CLOCK`
wall time, not `CPU_TIME`, and writes versioned outer and Schur CSV schemas.

The Schur rows record independent `reached_tolerance`, `hit_maxiter`,
`breakdown`, and `nonfinite` fields. Outer residual/iteration fields are left
missing when the hook does not own the outer Krylov state; no zero is used as
a fabricated residual. K accounting has separate primal block solve,
matrix-free Schur action, full upper correction, and setup/rebuild stages.
Per-call timers are summed; cumulative timers use their final sample and are
never summed again. GPU synchronization is reported missing until the
backend supplies a synchronization hook.

Use `scripts/analysis/summarize_block_schur_probe.py` to validate schemas and
produce `INCOMPLETE` for missing physical measurements. It reports actual
outer progress, K-action and wall-time normalizations, `log10` reductions,
and workload signatures for peer matching.

## Acceptance profiles

`scripts/analysis/evaluate_solver_acceptance.py` provides two profiles:

- `diagnostic` can pass numerical correctness while keeping production
  readiness false;
- `production` additionally requires primal agreement and TES temperature
  and current parity. Missing metrics are `INCOMPLETE`, never PASS.

Nonfinite values and breakdown are hard failures. The normwise backward error
definition used by the exact oracle is
`||A*x-b||2 / (||A||F*||x||2 + ||b||2)`. A relative residual above `1e-11`
with absolute residual at or below `1e-14` is reported as a numerical-floor
warning, not silently discarded.

## Runtime evidence

The installed `ElmerSolver 26.1-devel` was tried against all four SIFs with a
30-second per-case bound. Each exited during setup with
`Hypre requested but not compiled with!`; these are capability results, not
solver convergence results. Consequently no lower/full CPU/GPU comparison is
claimed and production readiness remains false. See
`artifacts/hypre_phase20/runtime_status_20260905.json` and the ignored runtime
logs for the exact attempts.

The independent exact SuperLU Schur oracle was rerun on the existing one-rank
matrix/RHS/MUMPS artifacts. It records absolute residual
`1.1384928931822508e-16`, relative residual `3.180826027618793e-11`, normwise
backward error `6.844108235359397e-19`, constraint residual
`1.1589146152036896e-24`, and primal agreement `2.2119885798515585e-6`.
The numerical diagnostic passes; TES physical parity is still missing, so the
acceptance artifact remains `INCOMPLETE` for production.

The native source implementation was applied to a clean Elmer source worktree
and built successfully under WSL. A non-MPI CPU build with
`WITH_SuperLU=TRUE` linked `ElmerSolver`; an HYPRE-enabled MPI build with
`WITH_Hypre=TRUE`, `/usr/include/hypre`, and `libHYPRE.so` linked
`ElmerSolver_mpi`. The reviewable patch was regenerated from that source diff
and passes `git apply --check` against the clean baseline.

The stock Windows install is `CAPABILITY MISSING` because it lacks HYPRE. A
HYPRE-enabled native MPI binary was then built and executed against the
existing one-rank mesh/restart. Both CPU probes completed with `ALL DONE` in
under the 240-second bound:

| probe | outer rows | Schur rows | K actions | upper correction | wall seconds |
|---|---:|---:|---:|---:|---:|
| lower CPU | 30 | 30 | 991 | 0 | 36.664 |
| full CPU | 30 | 30 | 1021 | 30 | 40.227 |

Both probes recorded 30 Schur iterations/maxiter hits per row and no
breakdown/nonfinite flags. The native `BlockMatrixPrec` hook does not own the
outer Krylov state, so outer solver iteration/residual fields are explicitly
empty; the summaries are therefore `INCOMPLETE` for outer convergence, never
a fabricated PASS. GPU is `CAPABILITY MISSING`: the host exposes an NVIDIA
device, but the linked HYPRE library is CPU-only and no CUDA/HIP HYPRE build
was available.
