# Phase24 Gate 3: HYPRE steady state

- Status: **FAIL**
- Case: `case_phase24_g3_mesh_hybrid_fullconf_h12_hybfull_h1_87ab7d`
- Source project: `D:\Github\TES-Programs\Elmer-Projects\artifacts\phase24_conformal_hybrid\project.json`
- Mesh: `mesh_hybrid_fullconf_h12`; mortar: `False`
- Initial state: `T_0`, with no restart input
- HYPRE: `mumps`, max 4000 iterations, tolerance `1e-10`
- HYPRE GMRES dimension: `100`
- HYPRE reuse: `False`; preconditioner lagging: `None`
- BoomerAMG strong threshold override: `None`

## Criteria

| criterion | result |
|---|---|
| `normal_exit` | `True` |
| `independent_initial_temperature` | `True` |
| `finite_observables` | `True` |
| `mumps_reference_recorded` | `True` |
| `comsol_current_within_0p6_percent` | `False` |
| `hypre_linear_residual_recorded` | `False` |
| `elmer_nonlinear_residual_recorded` | `True` |
| `tes_circuit_residual_recorded` | `True` |

- HYPRE current [µA]: `139.67634818614016`
- Current evidence source: `iteration_series` (a non-series source is not a converged steady result)
- Stage 11 MUMPS delta [%]: `0.00034192792670214965`
- COMSOL delta [%]: `2.3618186407806108`
- Stage 11 MUMPS reference [µA]: `139.67587059533173`
- COMSOL reference [µA]: `143.055049`
- HYPRE final linear residual: `None`
- Elmer nonlinear residual: `6.7274933e-09`
- TES circuit residual [W]: `-2.3689197345461714e-15`

Project: `D:\Github\TES-Programs\Elmer-Projects\artifacts\phase24_gate3_hypre_steady_hybfull_h12_mumps\phase24_gate3_hypre_steady.json`
Solver log: `D:\Github\TES-Programs\Elmer-Projects\results\case_phase24_g3_mesh_hybrid_fullconf_h12_hybfull_h1_87ab7d\solver.log`
Launcher log: `D:\Github\TES-Programs\Elmer-Projects\results\case_phase24_g3_mesh_hybrid_fullconf_h12_hybfull_h1_87ab7d\gate3_launcher.log`

Gate 3 failure is a stop condition; Gate 4 and later must not be started from this bundle.
