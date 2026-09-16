# Phase24 difference evidence

Generated (UTC): `2026-09-15T15:39:38.942349+00:00`
Git revision: `3540f7dd4a021ba5c60917e5d7dc7e7cd3fa7d57`

- Candidate: `BDF2 / HYPRE reuse ON` — `D:\Github\TES-Programs\Elmer-Projects\results\case_phase24_short_bdf2_reuse_on_1us\case_phase24_short_bdf2_reuse_on_1us_series.csv`
- Reference: `CPU/MUMPS Phase23` — `D:\Github\TES-Programs\Elmer-Projects\results\case_p19_pulse_phase23_tight\case_p19_pulse_phase23_tight_series.csv`
- Event: `0.02002 s`
- Common post-event window: `0..0.901 µs` (20 rows)
- Baselines: candidate `147.801439420 µA`, reference `147.375245042 µA`, offset `+0.426194378 µA`

## Difference onset

- First raw absolute difference ≥ `0.1 µA`: **0.011 µs**
- First baseline-corrected absolute difference ≥ `0.1 µA`: **3.63797881e-12 µs**
- First raw relative difference ≥ `0.1%`: **0.601 µs**
- Maximum raw difference: `-0.682342660 µA` at `0.901 µs`
- Maximum baseline-corrected difference: `-1.108537038 µA` at `0.901 µs`

## Checkpoints

| t from pulse (µs) | candidate (µA) | reference (µA) | raw Δ (µA) | corrected Δ (µA) |
|---:|---:|---:|---:|---:|
| 0.1 | 147.496189303 | 147.389028372 | +0.107160931 | -0.319033447 |
| 0.2 | 147.475445111 | 147.386253795 | +0.089191316 | -0.337003062 |
| 0.3 | 147.425217760 | 147.380147953 | +0.045069808 | -0.381124571 |
| 0.4 | 147.341609721 | 147.370795763 | -0.029186042 | -0.455380420 |
| 0.5 | 147.227404131 | 147.358717780 | -0.131313649 | -0.557508027 |
| 0.6 | 147.090694970 | 147.344626368 | -0.253931398 | -0.680125776 |
| 0.7 | 146.939239590 | 147.329182491 | -0.389942901 | -0.816137279 |
| 0.8 | 146.781606123 | 147.312935484 | -0.531329361 | -0.957523739 |
| 0.9 | 146.615465644 | 147.296312315 | -0.680846671 | -1.107041049 |

The raw difference includes the pre-pulse operating-point offset. The corrected difference subtracts each run's own mean over the configured pre-pulse baseline window.
The comparison is based on recorded output samples; it does not treat internal nonlinear/Krylov epochs as physical time steps.

Files: `summary.json`, `aligned_difference.csv`, `checkpoints.csv`, `comparison.png`
