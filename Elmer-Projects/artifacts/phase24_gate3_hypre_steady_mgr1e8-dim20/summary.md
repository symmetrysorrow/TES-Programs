# Phase24 Gate 3: HYPRE steady state

- Status: **FAIL**
- Case: `case_phase24_gate3_hypre_steady_prod_v2_mgr1e8-dim20`
- Mesh: `mesh_singlepixel_prod_v2`; mortar: `True`
- Initial state: `T_0`, with no restart input
- HYPRE: `iterative_hypre_flexgmres_mgr`, max 10000 iterations, tolerance `1e-08`
- HYPRE GMRES dimension: `20`
- HYPRE reuse: `False`; preconditioner lagging: `None`
- BoomerAMG strong threshold override: `None`

## Criteria

| criterion | result |
|---|---|
| `normal_exit` | `False` |
| `independent_initial_temperature` | `True` |
| `finite_observables` | `True` |
| `mumps_current_within_0p05_percent` | `False` |
| `comsol_current_within_0p6_percent` | `True` |
| `hypre_linear_residual_recorded` | `False` |
| `elmer_nonlinear_residual_recorded` | `False` |
| `tes_circuit_residual_recorded` | `True` |

- HYPRE current [µA]: `143.53734493231093`
- Current evidence source: `iteration_series` (a non-series source is not a converged steady result)
- Stage 11 MUMPS delta [%]: `0.16721019605360993`
- COMSOL delta [%]: `0.3371400979429485`
- HYPRE final linear residual: `None`
- Elmer nonlinear residual: `None`
- TES circuit residual [W]: `0.0`

Project: `D:\Github\TES-Programs\Elmer-Projects\artifacts\phase24_gate3_hypre_steady_mgr1e8-dim20\phase24_gate3_hypre_steady.json`
Solver log: `D:\Github\TES-Programs\Elmer-Projects\results\case_phase24_gate3_hypre_steady_prod_v2_mgr1e8-dim20\solver.log`
Launcher log: `D:\Github\TES-Programs\Elmer-Projects\results\case_phase24_gate3_hypre_steady_prod_v2_mgr1e8-dim20\gate3_launcher.log`

Gate 3 failure is a stop condition; Gate 4 and later must not be started from this bundle.
