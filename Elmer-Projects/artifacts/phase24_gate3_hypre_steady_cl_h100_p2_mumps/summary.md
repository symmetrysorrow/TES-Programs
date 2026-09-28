# Phase24 Gate 3: HYPRE steady state

- Status: **FAIL**
- Case: `case_phase24_g3_mesh_hybrid_comsollike_h100_p2_cl_h100_p2_b097b5`
- Source project: `D:\Github\TES-Programs\Elmer-Projects\artifacts\phase24_conformal_hybrid\project.json`
- Mesh: `mesh_hybrid_comsollike_h100_p2`; mortar: `False`
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

- HYPRE current [µA]: `142.0978843413425`
- Current evidence source: `iteration_series` (a non-series source is not a converged steady result)
- Stage 11 MUMPS delta [%]: `1.1669520011889305`
- COMSOL delta [%]: `0.6690883442062153`
- Stage 11 MUMPS reference [µA]: `143.775677487`
- COMSOL reference [µA]: `143.055049`
- HYPRE final linear residual: `None`
- Elmer nonlinear residual: `4.8600575e-09`
- TES circuit residual [W]: `-1.2384598521229088e-14`

Project: `D:\Github\TES-Programs\Elmer-Projects\artifacts\phase24_gate3_hypre_steady_cl_h100_p2_mumps\phase24_gate3_hypre_steady.json`
Solver log: `D:\Github\TES-Programs\Elmer-Projects\results\case_phase24_g3_mesh_hybrid_comsollike_h100_p2_cl_h100_p2_b097b5\solver.log`
Launcher log: `D:\Github\TES-Programs\Elmer-Projects\results\case_phase24_g3_mesh_hybrid_comsollike_h100_p2_cl_h100_p2_b097b5\gate3_launcher.log`

Gate 3 failure is a stop condition; Gate 4 and later must not be started from this bundle.
