# Phase24 same-path direct-MUMPS isolation

Generated: `2026-09-15T17:46:52.413817+00:00`

## Controlled change

- Mesh, restart, pulse, inner circuit, BDF2, timestep grid, nonlinear controls: unchanged from Phase24 HYPRE case
- Phase24 assembly/operator settings: unchanged
- Linear backend: `iterative_hypre_flexgmres_boomeramg` -> `mumps`

## Run result

- Process exit code: `0`
- `ALL DONE`: `True`
- Abort/allocation markers: `0`
- MUMPS markers: `1`
- Relative-change samples: `46`; zero values: `0`

The two comparison bundles separate backend mismatch from the older CPU/MUMPS reference's broader Phase23 differences.

Direct MUMPS vs CPU/MUMPS: `D:\Github\TES-Programs\Elmer-Projects\artifacts\comparison\phase24_same_path_mumps_vs_cpu_1us`
Direct MUMPS vs native HYPRE: `D:\Github\TES-Programs\Elmer-Projects\artifacts\comparison\phase24_same_path_mumps_vs_hypre_1us`
Solver log: `D:\Github\TES-Programs\Elmer-Projects\results\case_phase24_short_same_path_mumps_1us\solver.log`
Generated project: `D:\Github\TES-Programs\Elmer-Projects\artifacts\phase24_same_path_mumps_1us\phase24_same_path_mumps.json`
