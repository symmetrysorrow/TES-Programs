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
The historical `runtime_status_20260905.json` records the pre-CUDA-HYPRE
capability result. The current CUDA-HYPRE evidence is in
`hypre_cuda_build_evidence_20260905.txt`,
`elmer_phase20_gpu_link_evidence_20260905.txt`, and the lower/full Nsight
reports `nsys_phase20_lower_gpu_migrated2_20260905.nsys-rep` and
`nsys_phase20_full_gpu_migrated_20260905.nsys-rep`.

Those reports contain real CUDA kernel launches and memory transfers, and both
Phase20 lower/full runs reached `ALL DONE`. Both nevertheless reported outer
linear non-convergence, while the independent P19 smoke ended in
`MPI_ABORT`; therefore correctness is `FAIL` and performance readiness remains
`NO`. The outer Krylov state is not owned by the instrumentation hook, so
missing outer fields remain missing rather than being fabricated.
