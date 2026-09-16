# Phase24 HYPRE linear-tolerance diagnostic

Generated: `2026-09-16T03:53:29.711329+00:00`

## Fixed conditions

- Same mesh, restart, pulse, inner circuit, Phase24 assembly, BDF2, nonlinear controls, HYPRE reuse, and adaptive preconditioner lagging
- Same approximately `1 us` post-pulse window
- Only `Linear System Convergence Tolerance` changes

## Results

| condition | exit | ALL DONE | HYPRE solves | zero-iteration solves | max iterations |
|---|---:|---:|---:|---:|---:|
| HYPRE tolerance 5e-7 | 0 | True | 46 | 23 | 1258 |
| HYPRE tolerance 1e-10 | 1 | False | 0 | 0 | None |

5e-7 vs CPU/MUMPS: `not run`
1e-10 vs CPU/MUMPS: `not run`
1e-10 vs 5e-7: `not run`

Interpretation: if the strict tolerance moves the absolute current toward the direct/COMSOL value and changes the baseline-corrected waveform, linear convergence was insufficient. If it does not, investigate HYPRE vector transfer/scaling or the solver backend implementation.

Generated project: `D:\Github\TES-Programs\Elmer-Projects\artifacts\phase24_hypre_tolerance_diagnostic_1us\phase24_hypre_tolerance_diagnostic.json`
