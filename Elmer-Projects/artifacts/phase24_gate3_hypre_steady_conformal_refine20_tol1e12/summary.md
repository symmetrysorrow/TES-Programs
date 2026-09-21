# Phase24 Gate 3: HYPRE steady state

- Status: **FAIL**
- Case: `case_phase24_g3_conformal_refine20_conformal_refine20`
- Source project: `D:\Github\TES-Programs\Elmer-Projects\artifacts\comparison\nomortar_refinement_probe\project.json`
- Mesh: `mesh_singlepixel_conformal_gpu_refine20`; mortar: `False`
- Initial state: `T_0`, with no restart input
- HYPRE: `iterative_hypre_flexgmres_boomeramg`, max 4000 iterations, tolerance `1e-12`
- HYPRE GMRES dimension: `100`
- HYPRE reuse: `False`; preconditioner lagging: `None`
- BoomerAMG strong threshold override: `0.5`

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
- Stage 11 MUMPS delta [%]: `0.810883833092967`
- COMSOL delta [%]: `0.3371400979429485`
- Stage 11 MUMPS reference [µA]: `144.71078126230952`
- COMSOL reference [µA]: `143.055049`
- HYPRE final linear residual: `None`
- Elmer nonlinear residual: `None`
- TES circuit residual [W]: `0.0`

Project: `D:\Github\TES-Programs\Elmer-Projects\artifacts\phase24_gate3_hypre_steady_conformal_refine20_tol1e12\phase24_gate3_hypre_steady.json`
Solver log: `D:\Github\TES-Programs\Elmer-Projects\results\case_phase24_g3_conformal_refine20_conformal_refine20\solver.log`
Launcher log: `D:\Github\TES-Programs\Elmer-Projects\results\case_phase24_g3_conformal_refine20_conformal_refine20\gate3_launcher.log`

Gate 3 failure is a stop condition; Gate 4 and later must not be started from this bundle.
