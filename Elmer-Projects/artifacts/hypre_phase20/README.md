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
`MPI_ABORT`. The isolated no-mortar smoke subsequently passed its HYPRE linear
solve on both CPU and GPU, with solution relative L2 difference
`2.9639969303996504e-08`; its process exit code is still `1` because the
existing nonlinear/linear-system-abort status is emitted after the linear solve.
The no-mortar GPU Nsight report records 35,290 kernel launches. Phase20 lower
correctness remains `FAIL` because CPU and GPU share the same outer stagnation
pattern, and performance readiness remains `NO`. The outer Krylov state is not
owned by the instrumentation hook; the captured HUTI history is recorded in
`gpu_correctness_20260905/phase20_lower_outer_history.json` rather than being
fabricated in the hook.

The CPU-only Phase20 lower contract run adds
`cpu_lower_contract_baseline_schur_trace.csv`, which records 960 rows across 30
Schur solves: 30 initial rows, 900 Arnoldi estimates, and 30 true-residual
rows. Every solve has an explicit RHS-scaled threshold and
`stopping_reason=maxiter`; the trace summary is
`cpu_lower_contract_baseline_summary.json`. The outer/inner correspondence is
in `cpu_lower_outer_inner_alignment_20260905.json`.

The only tuning experiment in this continuation changed CPU lower
`Block Schur Max Iterations` from 30 to 60, leaving restart unchanged. It did
not produce outer convergence: all 30 solves still hit the cap. The comparison
is `cpu_lower_tuning_comparison_20260905.json`; no GPU retest or performance
measurement was started.
