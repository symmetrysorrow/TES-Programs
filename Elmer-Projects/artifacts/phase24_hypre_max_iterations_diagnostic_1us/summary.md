# Phase24 HYPRE maximum-iteration diagnostic

Generated: `2026-09-16T07:30:49.776075+00:00`

- Same Phase24/BDF2/HYPRE policy as the 1e-8 diagnostic
- Only `Linear System Max Iterations` changed to `10000`; tolerance remains `1e-8`
- Exit code: `1`
- ALL DONE: `False`
- Successful solve records: `0`
- Failure telemetry records: `2`

## Decision

The solver still reaches the 10000-iteration limit without satisfying 1e-8; this is residual stagnation or an effective convergence barrier, not merely the previous 2000-iteration cap.

## Failure telemetry

```json
[
  {
    "phase": "backend_status",
    "method": 901,
    "solve_status": 256,
    "iterations": 10000,
    "final_relative_residual": 1.4565552596879816e-07,
    "requested_tolerance": 1e-08,
    "max_iterations": 10000,
    "matrix_epoch": 2,
    "preconditioner_epoch": 2,
    "case": 0,
    "hypre_error": 256,
    "hypre_global_error": 256
  },
  {
    "phase": "backend_status",
    "method": 901,
    "solve_status": 256,
    "iterations": 10000,
    "final_relative_residual": 1.2678624053027087e-07,
    "requested_tolerance": 1e-08,
    "max_iterations": 10000,
    "matrix_epoch": 2,
    "preconditioner_epoch": 2,
    "case": 3,
    "hypre_error": 256,
    "hypre_global_error": 256
  }
]
```

Solver log: `D:\Github\TES-Programs\Elmer-Projects\results\case_phase24_short_hypre_tol1e8_max10000_1us\solver.log`
Launcher log: `D:\Github\TES-Programs\Elmer-Projects\results\case_phase24_short_hypre_tol1e8_max10000_1us\max_iterations_diagnostic_launcher.log`

1e-8 vs CPU/MUMPS: `D:\Github\TES-Programs\Elmer-Projects\artifacts\comparison\phase24_hypre_tol1e8_max10000_vs_cpu_1us`
1e-8 vs HYPRE 5e-7: `D:\Github\TES-Programs\Elmer-Projects\artifacts\comparison\phase24_hypre_tol1e8_max10000_vs_tol5e7_1us`
