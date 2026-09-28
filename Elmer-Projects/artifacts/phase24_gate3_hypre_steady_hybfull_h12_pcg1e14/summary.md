# Phase24 Gate 3: HYPRE steady state

- Status: **FAIL**
- Case: `case_phase24_g3_mesh_hybrid_fullconf_h12_hybfull_h1_a3b87c`
- Source project: `D:\Github\TES-Programs\Elmer-Projects\artifacts\phase24_conformal_hybrid\project.json`
- Mesh: `mesh_hybrid_fullconf_h12`; mortar: `False`
- Initial state: `T_0`, with no restart input
- HYPRE: `iterative_hypre_pcg_boomeramg`, max 4000 iterations, tolerance `1e-14`
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
| `comsol_current_within_0p6_percent` | `False` |
| `hypre_linear_residual_recorded` | `True` |
| `elmer_nonlinear_residual_recorded` | `True` |
| `tes_circuit_residual_recorded` | `True` |

- HYPRE current [µA]: `139.67529924436258`
- Current evidence source: `iteration_series` (a non-series source is not a converged steady result)
- Stage 11 MUMPS delta [%]: `0.0007509802419698803`
- COMSOL delta [%]: `2.362551884231237`
- Stage 11 MUMPS reference [µA]: `139.67634818614016`
- COMSOL reference [µA]: `143.055049`
- HYPRE final linear residual: `-1.0`
- Elmer nonlinear residual: `2.9034636e-09`
- TES circuit residual [W]: `-4.616507177290305e-15`

Project: `D:\Github\TES-Programs\Elmer-Projects\artifacts\phase24_gate3_hypre_steady_hybfull_h12_pcg1e14\phase24_gate3_hypre_steady.json`
Solver log: `D:\Github\TES-Programs\Elmer-Projects\results\case_phase24_g3_mesh_hybrid_fullconf_h12_hybfull_h1_a3b87c\solver.log`
Launcher log: `D:\Github\TES-Programs\Elmer-Projects\results\case_phase24_g3_mesh_hybrid_fullconf_h12_hybfull_h1_a3b87c\gate3_launcher.log`

Gate 3 failure is a stop condition; Gate 4 and later must not be started from this bundle.
