# Phase24 Gate 3: HYPRE steady state

- Status: **FAIL**
- Case: `case_phase24_gate3_hypre_steady_prod_v2_mgr-reuse`
- Mesh: `mesh_singlepixel_prod_v2`; mortar: `True`
- Initial state: `T_0`, with no restart input
- HYPRE: `iterative_hypre_flexgmres_mgr`, max 4000 iterations, tolerance `5e-07`
- HYPRE reuse: `True`; preconditioner lagging: `adaptive`

## Criteria

| criterion | result |
|---|---|
| `normal_exit` | `True` |
| `independent_initial_temperature` | `True` |
| `finite_observables` | `True` |
| `mumps_current_within_0p05_percent` | `False` |
| `comsol_current_within_0p6_percent` | `True` |
| `hypre_linear_residual_recorded` | `True` |
| `elmer_nonlinear_residual_recorded` | `True` |
| `tes_circuit_residual_recorded` | `True` |

- HYPRE current [µA]: `142.68588147365492`
- Current evidence source: `iteration_series` (a non-series source is not a converged steady result)
- Stage 11 MUMPS delta [%]: `0.7594182554532944`
- COMSOL delta [%]: `0.2580597671495541`
- HYPRE final linear residual: `4.78449e-07`
- Elmer nonlinear residual: `2.0`
- TES circuit residual [W]: `-8.831606285879533e-13`

Project: `D:\Github\TES-Programs\Elmer-Projects\artifacts\phase24_gate3_hypre_steady_mgr-reuse\phase24_gate3_hypre_steady.json`
Solver log: `D:\Github\TES-Programs\Elmer-Projects\results\case_phase24_gate3_hypre_steady_prod_v2_mgr-reuse\solver.log`
Launcher log: `D:\Github\TES-Programs\Elmer-Projects\results\case_phase24_gate3_hypre_steady_prod_v2_mgr-reuse\gate3_launcher.log`

Gate 3 failure is a stop condition; Gate 4 and later must not be started from this bundle.
