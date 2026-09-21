# Phase24 Gate 3: HYPRE steady state

- Status: **FAIL**
- Case: `case_phase24_gate3_hypre_steady_prod_v2_mgr10k`
- Mesh: `mesh_singlepixel_prod_v2`; mortar: `True`
- Initial state: `T_0`, with no restart input
- HYPRE: `iterative_hypre_flexgmres_mgr`, max 10000 iterations, tolerance `1e-08`

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

- HYPRE current [µA]: `147.82755084258986`
- Current evidence source: `iteration_series` (a non-series source is not a converged steady result)
- Stage 11 MUMPS delta [%]: `2.8167047221058636`
- COMSOL delta [%]: `3.336129606016116`
- HYPRE final linear residual: `9.99797e-09`
- Elmer nonlinear residual: `0.0014566465`
- TES circuit residual [W]: `3.0966277928615794e-12`

Project: `D:\Github\TES-Programs\Elmer-Projects\artifacts\phase24_gate3_hypre_steady_mgr10k\phase24_gate3_hypre_steady.json`
Solver log: `D:\Github\TES-Programs\Elmer-Projects\results\case_phase24_gate3_hypre_steady_prod_v2_mgr10k\solver.log`
Launcher log: `D:\Github\TES-Programs\Elmer-Projects\results\case_phase24_gate3_hypre_steady_prod_v2_mgr10k\gate3_launcher.log`

Gate 3 failure is a stop condition; Gate 4 and later must not be started from this bundle.
