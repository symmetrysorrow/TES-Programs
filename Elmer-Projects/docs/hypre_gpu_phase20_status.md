# Phase20 bounded block-Schur experiments

Updated: 2026-09-04

Phase19 matrix-free Schur correctness is treated as a separate result from
convergence and performance.  The same-binary K/Schur parity and short-Bt
finite-tail regression remain PASS.  The lower CPU production attempt remains
`INCOMPLETE` after about 2,700 seconds; its partial iterate is not a solution.
The most likely direct cause is algorithmic cost: the inner matrix-free Schur
GMRES repeatedly reaches its iteration limit while each action invokes a K
solve.  Faster K actions on a GPU cannot remove that nested-Krylov cost.

Lower and full are now peer candidates.  Phase20 always generates all four
bounded combinations:

| case | backend | factorization | outer limit | Schur tolerance / max / restart |
|---|---|---|---:|---|
| `case_p20_hypre_block_lower_cpu_probe` | CPU | lower | 15 | `1e-4 / 30 / 30` |
| `case_p20_hypre_block_lower_gpu_probe` | GPU | lower | 15 | `1e-4 / 30 / 30` |
| `case_p20_hypre_block_full_cpu_probe` | CPU | full | 15 | `1e-4 / 30 / 30` |
| `case_p20_hypre_block_full_gpu_probe` | GPU | full | 15 | `1e-4 / 30 / 30` |

All four use the same mesh, restart, first timestep, and basic tolerances.
The GPU cases are generated even on CPU-only hosts; a GPU runtime test is
environment-dependent and must be recorded as `SKIP`/`NOT RUN`, not treated as
a generator failure.

## Probe output

Set `Block Schur Probe = Logical True` only in a diagnostic case.  The Phase20
source instrumentation writes `<prefix>_outer.csv` and `<prefix>_schur.csv`.
The rows contain outer application count, Schur solve count, inner iterations,
initial/final residuals, tolerance/max-iteration status, K-action counts,
wall time, and accumulated K/action time.  Solver-reported outer residuals
remain separate because `BlockMatrixPrec` does not own the outer Krylov
residual; the Elmer log or outer solver artifact supplies that value.

Summarize the files with:

```text
python scripts/analysis/summarize_block_schur_probe.py \
  --outer-csv case_p20_hypre_block_lower_cpu_probe_outer.csv \
  --schur-csv case_p20_hypre_block_lower_cpu_probe_schur.csv \
  --output results/phase20_lower_cpu_probe.json
```

The summary reports residual reduction per outer step, per K action, and per
second.  The objective is total cost for useful outer progress, not the
smallest possible inner residual.

## Parameter sweep

Generate the representative 16-case matrix (lower/full × CPU/GPU × four
points), then render it with the normal case builder:

```text
python scripts/prep/prepare_phase20_schur_sweep.py
python sync_elmer_parameters.py elmer_project_hypre_gpu_phase20_sweep.json
```

The points are `(tol, max iterations, K AMG cycles)` =
`(1e-2,5,1)`, `(1e-3,10,1)`, `(1e-4,20,1)`, `(1e-4,30,2)`.  Add `--all` for
the complete 3×4×2 matrix per factorization/backend.

## Acceptance

`scripts/analysis/evaluate_solver_acceptance.py` evaluates the metrics
separately.  Nonfinite values and breakdown are hard failures.  Production
numerical gates are absolute original-system residual (`1e-12`), backward
error (`1e-12`), constraint absolute residual (`1e-12`), and primal comparison
against MUMPS (`1e-4` when available).  Relative residual remains a reported
diagnostic.  If it exceeds `1e-11` while the absolute residual is at or below
the `1e-14` numerical-floor warning level, it is explicitly reported as a
floor warning rather than a sole hard failure.  Missing primary metrics are
`INCOMPLETE`; they are never inferred as zero.

Physical TES temperature/current parity is required for production promotion
once those measurements are available.  A passing Python test or generated
SIF proves configuration integrity only; it does not prove that an Elmer
solver converges.

## Explicit/reusable Schur prototype

`scripts/analysis/compare_schur_strategies.py` records setup and solve time,
Schur size, backward error, and the reuse boundary.  Its exact/ILU results
build the factorization once and reuse it for all columns of `B^T` and the
final primal action within one matrix/timestep.  The 2,898-DOF constraint
block and approximately 67 MB dense Schur observation justify evaluating a
frozen or sparse explicit Schur candidate.  Reuse across timesteps remains an
experiment, not an assumption: matrix changes from nonlinear/timestep
assembly must be checked before reusing a factorization.

## Current Phase20 runtime status

The four bounded SIFs and the representative sweep generator are implemented.
The modified source compiles and links.  A CPU lower probe was attempted with
a 120-second limit and reached Elmer startup, HYPRE setup, matrix dump, and
Schur diagnostic logging, but did not finish the large mounted-workspace dump;
it is recorded as `INCOMPLETE` with no convergence metrics and no accepted
iterate.  Full CPU and both GPU probes remain `NOT RUN`.  The next runtime
action is to avoid or relocate the large dump, then run all four probes under
the same one-rank environment and compare their JSON summaries.
