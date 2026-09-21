# Phase24 Gate 3: HYPRE steady state

- Status: **PASS**
- Case: `case_phase24_g3_s32m_s32m2_10f17b`
- Source project: `D:\Github\TES-Programs\Elmer-Projects\elmer_project_singlepixel_conformal_gpu.json`
- Mesh: `mesh_singlepixel_gpu_fine_stycast32_mortar`; mortar: `True`
- Initial state: `T_0`, with no restart input
- HYPRE: `iterative_hypre_flexgmres_boomeramg`, max 10000 iterations, tolerance `2e-06`
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
| `comsol_current_within_0p6_percent` | `True` |
| `hypre_linear_residual_recorded` | `True` |
| `elmer_nonlinear_residual_recorded` | `True` |
| `tes_circuit_residual_recorded` | `True` |

- HYPRE current [µA]: `143.56758901897666`
- Current evidence source: `iteration_series` (a non-series source is not a converged steady result)
- Stage 11 MUMPS delta [%]: `0.023103734455335555`
- COMSOL delta [%]: `0.3582816702797153`
- Stage 11 MUMPS reference [µA]: `143.534427206063`
- COMSOL reference [µA]: `143.055049`
- HYPRE final linear residual: `1.99717e-06`
- Elmer nonlinear residual: `3.9939422e-12`
- TES circuit residual [W]: `4.845735344315754e-14`

Project: `D:\Github\TES-Programs\Elmer-Projects\artifacts\phase24_gate3_hypre_steady_s32m2\phase24_gate3_hypre_steady.json`
Solver log: `D:\Github\TES-Programs\Elmer-Projects\results\case_phase24_g3_s32m_s32m2_10f17b\solver.log`
Launcher log: `D:\Github\TES-Programs\Elmer-Projects\results\case_phase24_g3_s32m_s32m2_10f17b\gate3_launcher.log`

Gate 3 passed with COMSOL as the primary physical target; the MUMPS delta is retained as a backend diagnostic. Gate 4 and later must use this no-mortar policy and remain subject to their own gates.
