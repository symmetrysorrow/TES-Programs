# Phase24 Gate 6: HIP HYPRE GPU parity

- Status: **PASS**
- Mesh: `mesh_singlepixel_conformal_gpu_refine20`; mortar: `False`
- Backend: `HIP HYPRE device`
- Steady GPU current [µA]: `143.52217883888775`
- CPU HYPRE steady current [µA]: `143.52217883888775`
- Steady delta [%]: `0.0`

| criterion | result |
|---|---|
| `steady_gpu_normal_exit` | `True` |
| `transient_gpu_normal_exit` | `True` |
| `steady_current_within_cpu_hypre_0p05_percent` | `True` |
| `transient_waveform_max_within_gate4_0p05_uA` | `True` |
| `transient_waveform_rmse_within_gate4_0p03_uA` | `True` |
| `gpu_execution_marker_present` | `True` |
| `no_cpu_fallback` | `True` |
| `no_device_error_or_native_crash` | `True` |
| `same_mesh_no_mortar_and_gpu_sif` | `True` |

- GPU solver log (steady): `D:\Github\TES-Programs\Elmer-Projects\results\case_phase24_g6_gpu_steady_conformal_refine20\solver.log`
- GPU solver log (transient): `D:\Github\TES-Programs\Elmer-Projects\results\case_phase24_g6_gpu_short_conformal_refine20\solver.log`
- Machine-readable artifact: `D:\Github\TES-Programs\Elmer-Projects\artifacts\phase24_gate6_gpu_parity\summary.json`
