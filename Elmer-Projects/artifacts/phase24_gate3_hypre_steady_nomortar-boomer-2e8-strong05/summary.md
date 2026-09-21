# Phase24 Gate 3: HYPRE steady state

- Status: **PASS**
- Case: `case_phase24_gate3_hypre_steady_prod_v2_nomortar-boomer-2e8-strong05`
- Mesh: `mesh_singlepixel_prod_v2`; mortar: `False`
- Initial state: `T_0`, with no restart input
- HYPRE: `iterative_hypre_flexgmres_boomeramg`, max 4000 iterations, tolerance `2e-08`
- HYPRE GMRES dimension: `100`
- HYPRE reuse: `False`; preconditioner lagging: `None`
- BoomerAMG strong threshold override: `0.5`

## Criteria

| criterion | result |
|---|---|
| `normal_exit` | `True` |
| `independent_initial_temperature` | `True` |
| `finite_observables` | `True` |
| `mumps_current_within_0p05_percent` | `True` |
| `comsol_current_within_0p6_percent` | `True` |
| `hypre_linear_residual_recorded` | `True` |
| `elmer_nonlinear_residual_recorded` | `True` |
| `tes_circuit_residual_recorded` | `True` |

- HYPRE current [µA]: `143.5373449311356`
- Current evidence source: `iteration_series` (a non-series source is not a converged steady result)
- Stage 11 MUMPS delta [%]: `8.188272889050257e-10`
- COMSOL delta [%]: `0.33714009712136067`
- HYPRE final linear residual: `1.97505e-08`
- Elmer nonlinear residual: `0.0`
- TES circuit residual [W]: `-1.961658822669505e-21`

Project: `D:\Github\TES-Programs\Elmer-Projects\artifacts\phase24_gate3_hypre_steady_nomortar-boomer-2e8-strong05\phase24_gate3_hypre_steady.json`
Solver log: `D:\Github\TES-Programs\Elmer-Projects\results\case_phase24_gate3_hypre_steady_prod_v2_nomortar-boomer-2e8-strong05\solver.log`
Launcher log: `D:\Github\TES-Programs\Elmer-Projects\results\case_phase24_gate3_hypre_steady_prod_v2_nomortar-boomer-2e8-strong05\gate3_launcher.log`

Gate 3 passed; Gate 4 and later must use this no-mortar policy and remain subject to their own gates.
