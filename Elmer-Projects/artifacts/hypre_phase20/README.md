# Phase20 probe artifacts

This directory holds machine-readable summaries from bounded lower/full
CPU/GPU probes.  A probe is intentionally limited to 15 outer iterations and
is not a production acceptance run.

Expected files are `<case>_outer.csv`, `<case>_schur.csv`, and the JSON summary
produced by `summarize_block_schur_probe.py`. GPU files are optional when no
CUDA/HIP runtime is available. Missing or unavailable runtime artifacts must
be labelled with an explicit status such as `CAPABILITY MISSING`, `NOT RUN`,
or `SKIP`; they must not be interpreted as solver convergence evidence.

The evaluator is `scripts/analysis/evaluate_solver_acceptance.py`.  It keeps
absolute residual, backward error, constraint residual, relative residual,
primal comparison, and physical TES parity as separate fields.

`exact_schur_recomputed.json` and `exact_schur_acceptance.json` contain the
rerun SuperLU oracle; its backward error is a real normwise value, not null.
`runtime_status_20260905.json` records the four bounded attempts. The
installed Elmer binary lacks HYPRE, so those attempts are `CAPABILITY MISSING`
for solver comparison purposes and do not provide convergence evidence. A
separate native source worktree was built with both CPU/SuperLU and HYPRE/MPI
configurations. The native HYPRE binary completed the bounded lower/full CPU
probes; GPU remains capability-missing and outer Krylov residuals are not
owned by the instrumentation hook, so the overall artifact stays `INCOMPLETE`.
