# Phase24 Gate 3: HYPRE steady state

- Status: **FAIL**
- Case: `case_phase24_g3_s32m_s32m_dea252`
- Source project: `D:\Github\TES-Programs\Elmer-Projects\elmer_project_singlepixel_conformal_gpu.json`
- Mesh: `mesh_singlepixel_gpu_fine_stycast32_mortar`; mortar: `True`
- Initial state: `T_0`, with no restart input
- HYPRE: `iterative_hypre_flexgmres_boomeramg`, max 4000 iterations, tolerance `2e-08`
- HYPRE GMRES dimension: `100`
- HYPRE reuse: `False`; preconditioner lagging: `None`
- BoomerAMG strong threshold override: `None`

## Criteria

| criterion | result |
|---|---|
| `normal_exit` | `False` |
| `independent_initial_temperature` | `True` |
| `finite_observables` | `True` |
| `mumps_reference_recorded` | `True` |
| `comsol_current_within_0p6_percent` | `True` |
| `hypre_linear_residual_recorded` | `True` |
| `elmer_nonlinear_residual_recorded` | `False` |
| `tes_circuit_residual_recorded` | `True` |

- HYPRE current [µA]: `143.53734493231093`
- Current evidence source: `iteration_series` (a non-series source is not a converged steady result)
- Stage 11 MUMPS delta [%]: `0.0020327710255566987`
- COMSOL delta [%]: `0.3371400979429485`
- Stage 11 MUMPS reference [µA]: `143.534427206063`
- COMSOL reference [µA]: `143.055049`
- HYPRE final linear residual: `1.8012078650071562e-06`
- Elmer nonlinear residual: `None`
- TES circuit residual [W]: `0.0`

Project: `D:\Github\TES-Programs\Elmer-Projects\artifacts\phase24_gate3_hypre_steady_s32m\phase24_gate3_hypre_steady.json`
Solver log: `D:\Github\TES-Programs\Elmer-Projects\results\case_phase24_g3_s32m_s32m_dea252\solver.log`
Launcher log: `D:\Github\TES-Programs\Elmer-Projects\results\case_phase24_g3_s32m_s32m_dea252\gate3_launcher.log`

Gate 3 failure is a stop condition; Gate 4 and later must not be started from this bundle.
