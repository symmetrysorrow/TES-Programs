# Phase24 HYPRE failure telemetry diagnostic

Generated: `2026-09-16T05:59:55.770220+00:00`

- Case: `case_phase24_short_hypre_tol1e8_1us`
- Numerical policy: unchanged; only failure telemetry and capture suppression are diagnostic changes
- Exit code: `1`
- ALL DONE: `False`
- Telemetry records: `2`
- Telemetry missing: `False`

## Native failure records

```json
[
  {
    "phase": "backend_status",
    "method": 901,
    "solve_status": 256,
    "iterations": 2000,
    "final_relative_residual": 1.8283706726485033e-07,
    "requested_tolerance": 1e-08,
    "max_iterations": 2000,
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
    "iterations": 2000,
    "final_relative_residual": 1.6395401610311612e-07,
    "requested_tolerance": 1e-08,
    "max_iterations": 2000,
    "matrix_epoch": 2,
    "preconditioner_epoch": 2,
    "case": 3,
    "hypre_error": 256,
    "hypre_global_error": 256
  }
]
```

Solver log: `D:\Github\TES-Programs\Elmer-Projects\results\case_phase24_short_hypre_tol1e8_1us\solver.log`
Launcher log: `D:\Github\TES-Programs\Elmer-Projects\results\case_phase24_short_hypre_tol1e8_1us\failure_telemetry_launcher.log`
