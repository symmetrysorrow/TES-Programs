# Phase24 difference evidence

Generated (UTC): `2026-09-15T17:33:29.630859+00:00`
Git revision: `45778dd1f918b861eb6b7b342672336f7c28b5bd`

- Candidate: `native HYPRE / reuse ON / lagging disabled` — `D:\Github\TES-Programs\Elmer-Projects\results\case_phase24_short_bdf2_reuse_on_lagging_disabled_1us\case_phase24_short_bdf2_reuse_on_lagging_disabled_1us_series.csv`
- Reference: `CPU/MUMPS Phase23` — `D:\Github\TES-Programs\Elmer-Projects\results\case_p19_pulse_phase23_tight\case_p19_pulse_phase23_tight_series.csv`
- Event: `0.02002 s`
- Common post-event window: `0..0.901 µs` (20 rows)
- Baselines: candidate `147.787051862 µA`, reference `147.375245042 µA`, offset `+0.411806820 µA`

## Difference onset

- First raw absolute difference ≥ `0.1 µA`: **3.63797881e-12 µs**
- First baseline-corrected absolute difference ≥ `0.1 µA`: **3.63797881e-12 µs**
- First raw relative difference ≥ `0.1%`: **3.63797881e-12 µs**
- Maximum raw difference: `-0.734105542 µA` at `0.901 µs`
- Maximum baseline-corrected difference: `-1.145912362 µA` at `0.901 µs`

## Checkpoints

| t from pulse (µs) | candidate (µA) | reference (µA) | raw Δ (µA) | corrected Δ (µA) |
|---:|---:|---:|---:|---:|
| 0.1 | 147.541519940 | 147.389028372 | +0.152491568 | -0.259315252 |
| 0.2 | 147.472097385 | 147.386253795 | +0.085843591 | -0.325963229 |
| 0.3 | 147.394469311 | 147.380147953 | +0.014321359 | -0.397485461 |
| 0.4 | 147.296508109 | 147.370795763 | -0.074287654 | -0.486094474 |
| 0.5 | 147.177038706 | 147.358717780 | -0.181679074 | -0.593485894 |
| 0.6 | 147.038104823 | 147.344626368 | -0.306521545 | -0.718328365 |
| 0.7 | 146.886196701 | 147.329182491 | -0.442985790 | -0.854792610 |
| 0.8 | 146.730576525 | 147.312935484 | -0.582358958 | -0.994165778 |
| 0.9 | 146.563710299 | 147.296312315 | -0.732602016 | -1.144408836 |

The raw difference includes the pre-pulse operating-point offset. The corrected difference subtracts each run's own mean over the configured pre-pulse baseline window.
The comparison is based on recorded output samples; it does not treat internal nonlinear/Krylov epochs as physical time steps.

Files: `summary.json`, `aligned_difference.csv`, `checkpoints.csv`, `comparison.png`
