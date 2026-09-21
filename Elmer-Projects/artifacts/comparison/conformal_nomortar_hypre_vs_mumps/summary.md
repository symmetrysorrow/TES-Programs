# Phase24 difference evidence

Generated (UTC): `2026-09-17T04:56:09.871241+00:00`
Git revision: `0081cc2e7a866020b7931bfb488259e9c02b9f29`

- Candidate: `conformal no-mortar native HYPRE` — `D:\Github\TES-Programs\Elmer-Projects\results\case_tes_pulse_singlepixel_conformal_gpu_fine_hypre_nomortar_smoke_5stage\case_tes_pulse_singlepixel_conformal_gpu_fine_hypre_nomortar_smoke_5stage_series.csv`
- Reference: `conformal no-mortar MUMPS` — `D:\Github\TES-Programs\Elmer-Projects\results\case_tes_pulse_singlepixel_conformal_gpu_fine_hybrid_177step\case_tes_pulse_singlepixel_conformal_gpu_fine_hybrid_177step_series.csv`
- Event: `0.02002 s`
- Common post-event window: `0..0.901 µs` (20 rows)
- Baselines: candidate `154.851462193 µA`, reference `154.851454831 µA`, offset `+0.000007362 µA`

## Difference onset

- First raw absolute difference ≥ `0.1 µA`: **n/a µs**
- First baseline-corrected absolute difference ≥ `0.1 µA`: **n/a µs**
- First raw relative difference ≥ `0.1%`: **n/a µs**
- Maximum raw difference: `+0.074747817 µA` at `0.901 µs`
- Maximum baseline-corrected difference: `+0.074740455 µA` at `0.901 µs`

## Checkpoints

| t from pulse (µs) | candidate (µA) | reference (µA) | raw Δ (µA) | corrected Δ (µA) |
|---:|---:|---:|---:|---:|
| 0.001 | 154.851435768 | 154.851352403 | +0.000083364 | +0.000076003 |
| 0.01 | 154.851435728 | 154.851349730 | +0.000085998 | +0.000078636 |
| 0.05 | 154.851433925 | 154.851334560 | +0.000099365 | +0.000092003 |
| 0.1 | 154.851415579 | 154.851216625 | +0.000198954 | +0.000191592 |
| 0.2 | 154.851369697 | 154.850881600 | +0.000488097 | +0.000480735 |
| 0.3 | 154.851204476 | 154.847602926 | +0.003601550 | +0.003594188 |
| 0.5 | 154.850797027 | 154.831393799 | +0.019403228 | +0.019395866 |
| 0.7 | 154.850519032 | 154.805781928 | +0.044737104 | +0.044729742 |
| 0.9 | 154.850400671 | 154.775805230 | +0.074595441 | +0.074588079 |

The raw difference includes the pre-pulse operating-point offset. The corrected difference subtracts each run's own mean over the configured pre-pulse baseline window.
The comparison is based on recorded output samples; it does not treat internal nonlinear/Krylov epochs as physical time steps.

Files: `summary.json`, `aligned_difference.csv`, `checkpoints.csv`, `comparison.png`
