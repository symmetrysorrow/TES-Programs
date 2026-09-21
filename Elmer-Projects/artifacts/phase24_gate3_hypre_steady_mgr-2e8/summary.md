# Phase24 Gate 3: HYPRE steady state

- Status: **FAIL**
- Case: `case_phase24_gate3_hypre_steady_prod_v2_mgr-2e8`
- Mesh: `mesh_singlepixel_prod_v2`; mortar: `True`
- Initial state: `T_0`, with no restart input
- HYPRE: `iterative_hypre_flexgmres_mgr`, max 4000 iterations, tolerance `2e-08`
- HYPRE reuse: `False`; preconditioner lagging: `None`
- BoomerAMG strong threshold override: `None`

## Criteria

| criterion | result |
|---|---|
| `normal_exit` | `False` |
| `independent_initial_temperature` | `True` |
| `finite_observables` | `True` |
| `mumps_current_within_0p05_percent` | `False` |
| `comsol_current_within_0p6_percent` | `False` |
| `hypre_linear_residual_recorded` | `True` |
| `elmer_nonlinear_residual_recorded` | `True` |
| `tes_circuit_residual_recorded` | `True` |

- HYPRE current [µA]: `140.96761167892606`
- Current evidence source: `iteration_series` (a non-series source is not a converged steady result)
- Stage 11 MUMPS delta [%]: `1.9545056198219903`
- COMSOL delta [%]: `1.459184653506312`
- HYPRE final linear residual: `1.98036e-08`
- Elmer nonlinear residual: `0.00013734156`
- TES circuit residual [W]: `-2.4774913271676197e-12`

Project: `D:\Github\TES-Programs\Elmer-Projects\artifacts\phase24_gate3_hypre_steady_mgr-2e8\phase24_gate3_hypre_steady.json`
Solver log: `D:\Github\TES-Programs\Elmer-Projects\results\case_phase24_gate3_hypre_steady_prod_v2_mgr-2e8\solver.log`
Launcher log: `D:\Github\TES-Programs\Elmer-Projects\results\case_phase24_gate3_hypre_steady_prod_v2_mgr-2e8\gate3_launcher.log`

Gate 3 failure is a stop condition; Gate 4 and later must not be started from this bundle.
