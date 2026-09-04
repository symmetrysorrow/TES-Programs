# Phase20 probe artifacts

This directory holds machine-readable summaries from bounded lower/full
CPU/GPU probes.  A probe is intentionally limited to 15 outer iterations and
is not a production acceptance run.

Expected files are `<case>_outer.csv`, `<case>_schur.csv`, and the JSON summary
produced by `summarize_block_schur_probe.py`.  GPU files are optional when no
CUDA/HIP runtime is available.  Missing runtime artifacts must be labelled
`NOT RUN` or `SKIP`.

The evaluator is `scripts/analysis/evaluate_solver_acceptance.py`.  It keeps
absolute residual, backward error, constraint residual, relative residual,
primal comparison, and physical TES parity as separate fields.

`exact_schur_recomputed.json` and `exact_schur_acceptance.json` contain the
rerun SuperLU oracle; its backward error is a real normwise value, not null.
`runtime_status_20260905.json` records the four bounded attempts.  The
installed Elmer binary lacks HYPRE, so those attempts are `NOT RUN` for solver
comparison purposes and do not provide convergence evidence.
