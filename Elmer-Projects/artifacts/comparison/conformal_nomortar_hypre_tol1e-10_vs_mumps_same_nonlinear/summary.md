# Phase24 difference evidence

Generated (UTC): `2026-09-17T05:20:55.139420+00:00`
Git revision: `b98fb4bce552010dec9f3c8a771fb64ec2e6c83c`

- Candidate: `conformal no-mortar native HYPRE tol=1e-10` — `D:\Github\TES-Programs\Elmer-Projects\results\case_tes_pulse_singlepixel_conformal_gpu_fine_hypre_nomortar_smoke_5stage_tol1e-10\case_tes_pulse_singlepixel_conformal_gpu_fine_hypre_nomortar_smoke_5stage_tol1e-10_series.csv`
- Reference: `conformal no-mortar native MUMPS same nonlinear tolerance` — `D:\Github\TES-Programs\Elmer-Projects\results\case_tes_pulse_singlepixel_conformal_gpu_fine_mumps_nomortar_smoke_5stage_tol1e-10\case_tes_pulse_singlepixel_conformal_gpu_fine_mumps_nomortar_smoke_5stage_tol1e-10_series.csv`
- Event: `0.02002 s`
- Common post-event window: `0..0.901 µs` (20 rows)
- Baselines: candidate `154.851462193 µA`, reference `154.851285427 µA`, offset `+0.000176766 µA`

## Difference onset

- First raw absolute difference ≥ `0.1 µA`: **n/a µs**
- First baseline-corrected absolute difference ≥ `0.1 µA`: **n/a µs**
- First raw relative difference ≥ `0.1%`: **n/a µs**
- Maximum raw difference: `+0.013549053 µA` at `0.901 µs`
- Maximum baseline-corrected difference: `+0.013372287 µA` at `0.901 µs`

## Checkpoints

| t from pulse (µs) | candidate (µA) | reference (µA) | raw Δ (µA) | corrected Δ (µA) |
|---:|---:|---:|---:|---:|
| 0.001 | 154.851435768 | 154.850923810 | +0.000511958 | +0.000335192 |
| 0.01 | 154.851435728 | 154.850917935 | +0.000517793 | +0.000341027 |
| 0.05 | 154.851433910 | 154.850888749 | +0.000545161 | +0.000368395 |
| 0.1 | 154.851415261 | 154.850760951 | +0.000654310 | +0.000477544 |
| 0.2 | 154.851368524 | 154.850417011 | +0.000951513 | +0.000774748 |
| 0.3 | 154.849405227 | 154.847133342 | +0.002271886 | +0.002095120 |
| 0.5 | 154.836736681 | 154.830919600 | +0.005817081 | +0.005640315 |
| 0.7 | 154.815078369 | 154.805306012 | +0.009772356 | +0.009595591 |
| 0.9 | 154.788859731 | 154.775328454 | +0.013531278 | +0.013354512 |

The raw difference includes the pre-pulse operating-point offset. The corrected difference subtracts each run's own mean over the configured pre-pulse baseline window.
The comparison is based on recorded output samples; it does not treat internal nonlinear/Krylov epochs as physical time steps.

Files: `summary.json`, `aligned_difference.csv`, `checkpoints.csv`, `comparison.png`
