# Phase24 difference evidence

Generated (UTC): `2026-09-17T05:09:15.791822+00:00`
Git revision: `b98fb4bce552010dec9f3c8a771fb64ec2e6c83c`

- Candidate: `conformal no-mortar native HYPRE tol=1e-10` — `D:\Github\TES-Programs\Elmer-Projects\results\case_tes_pulse_singlepixel_conformal_gpu_fine_hypre_nomortar_smoke_5stage_tol1e-10\case_tes_pulse_singlepixel_conformal_gpu_fine_hypre_nomortar_smoke_5stage_tol1e-10_series.csv`
- Reference: `conformal no-mortar MUMPS` — `D:\Github\TES-Programs\Elmer-Projects\results\case_tes_pulse_singlepixel_conformal_gpu_fine_hybrid_177step\case_tes_pulse_singlepixel_conformal_gpu_fine_hybrid_177step_series.csv`
- Event: `0.02002 s`
- Common post-event window: `0..0.901 µs` (20 rows)
- Baselines: candidate `154.851462193 µA`, reference `154.851454831 µA`, offset `+0.000007362 µA`

## Difference onset

- First raw absolute difference ≥ `0.1 µA`: **n/a µs**
- First baseline-corrected absolute difference ≥ `0.1 µA`: **n/a µs**
- First raw relative difference ≥ `0.1%`: **n/a µs**
- Maximum raw difference: `+0.013072272 µA` at `0.901 µs`
- Maximum baseline-corrected difference: `+0.013064911 µA` at `0.901 µs`

## Checkpoints

| t from pulse (µs) | candidate (µA) | reference (µA) | raw Δ (µA) | corrected Δ (µA) |
|---:|---:|---:|---:|---:|
| 0.001 | 154.851435768 | 154.851352403 | +0.000083364 | +0.000076003 |
| 0.01 | 154.851435728 | 154.851349730 | +0.000085998 | +0.000078636 |
| 0.05 | 154.851433910 | 154.851334560 | +0.000099350 | +0.000091988 |
| 0.1 | 154.851415261 | 154.851216625 | +0.000198636 | +0.000191274 |
| 0.2 | 154.851368524 | 154.850881600 | +0.000486925 | +0.000479563 |
| 0.3 | 154.849405227 | 154.847602926 | +0.001802301 | +0.001794939 |
| 0.5 | 154.836736681 | 154.831393799 | +0.005342883 | +0.005335521 |
| 0.7 | 154.815078369 | 154.805781928 | +0.009296441 | +0.009289079 |
| 0.9 | 154.788859731 | 154.775805230 | +0.013054501 | +0.013047139 |

The raw difference includes the pre-pulse operating-point offset. The corrected difference subtracts each run's own mean over the configured pre-pulse baseline window.
The comparison is based on recorded output samples; it does not treat internal nonlinear/Krylov epochs as physical time steps.

Files: `summary.json`, `aligned_difference.csv`, `checkpoints.csv`, `comparison.png`
