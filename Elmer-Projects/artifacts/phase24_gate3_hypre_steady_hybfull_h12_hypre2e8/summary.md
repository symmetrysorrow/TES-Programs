# Phase24 Gate 3: HYPRE steady state

- Status: **PASS**
- Case: `case_phase24_g3_mesh_hybrid_fullconf_h12_hybfull_h1_e8b223`
- Source project: `D:\Github\TES-Programs\Elmer-Projects\artifacts\phase24_conformal_hybrid\project.json`
- Mesh: `mesh_hybrid_fullconf_h12`; mortar: `False`
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
| `mumps_reference_recorded` | `True` |
| `comsol_current_within_0p6_percent` | `True` |
| `hypre_linear_residual_recorded` | `True` |
| `elmer_nonlinear_residual_recorded` | `True` |
| `tes_circuit_residual_recorded` | `True` |

- HYPRE current [µA]: `142.9612095288064`
- Current evidence source: `iteration_series` (a non-series source is not a converged steady result)
- Stage 11 MUMPS delta [%]: `2.3521163100482334`
- COMSOL delta [%]: `0.06559675582901246`
- Stage 11 MUMPS reference [µA]: `139.67587059533173`
- COMSOL reference [µA]: `143.055049`
- HYPRE final linear residual: `1.88722e-08`
- Elmer nonlinear residual: `0.0`
- TES circuit residual [W]: `-9.262924380105605e-13`

Project: `D:\Github\TES-Programs\Elmer-Projects\artifacts\phase24_gate3_hypre_steady_hybfull_h12_hypre2e8\phase24_gate3_hypre_steady.json`
Solver log: `D:\Github\TES-Programs\Elmer-Projects\results\case_phase24_g3_mesh_hybrid_fullconf_h12_hybfull_h1_e8b223\solver.log`
Launcher log: `D:\Github\TES-Programs\Elmer-Projects\results\case_phase24_g3_mesh_hybrid_fullconf_h12_hybfull_h1_e8b223\gate3_launcher.log`

Gate 3 passed with COMSOL as the primary physical target; the MUMPS delta is retained as a backend diagnostic. Gate 4 and later must use this no-mortar policy and remain subject to their own gates.
