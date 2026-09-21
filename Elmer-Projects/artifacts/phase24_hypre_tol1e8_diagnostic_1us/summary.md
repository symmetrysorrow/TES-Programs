# Phase24 HYPRE intermediate tolerance diagnostic

Generated: `2026-09-16T05:16:54.691463+00:00`

## Fixed conditions

- Same mesh, restart, pulse, inner circuit, Phase24 assembly, BDF2, nonlinear controls, HYPRE reuse, and adaptive preconditioner lagging
- Same approximately `1 us` post-pulse window
- Only `Linear System Convergence Tolerance` changes from the previous `5e-7` case

## Result

| condition | exit | ALL DONE | HYPRE solves | zero-iteration solves | max iterations | tolerance failure |
|---|---:|---:|---:|---:|---:|---:|
| HYPRE tolerance 1e-8 | 1 | False | 0 | 0 | None | True |

## Comparisons

- 1e-8 vs CPU/MUMPS: `D:\Github\TES-Programs\Elmer-Projects\artifacts\comparison\phase24_hypre_tol1e8_vs_cpu_1us`
- 1e-8 vs HYPRE 5e-7: `D:\Github\TES-Programs\Elmer-Projects\artifacts\comparison\phase24_hypre_tol1e8_vs_tol5e7_1us`
- Existing HYPRE 5e-7 vs CPU/MUMPS: `D:\Github\TES-Programs\Elmer-Projects\artifacts\comparison\phase24_hypre_tol5e7_vs_cpu_1us`

Interpretation: compare both absolute current and baseline-corrected waveform. A successful 1e-8 run that moves both toward the reference supports insufficient linear convergence as a cause; unchanged results point to another backend or scaling issue. The selector-mismatch warning, if present, is recorded separately and is not treated as the root cause.

Generated project: `D:\Github\TES-Programs\Elmer-Projects\artifacts\phase24_hypre_tol1e8_diagnostic_1us\phase24_hypre_tol1e8_diagnostic.json`
