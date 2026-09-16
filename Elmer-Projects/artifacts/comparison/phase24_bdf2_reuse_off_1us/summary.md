# Phase24 difference evidence

Generated (UTC): `2026-09-15T16:13:40.839399+00:00`
Git revision: `d3899fbeaf931f89afdf154fa43ca403cc15e61d`

- Candidate: `BDF2 / reuse OFF` — `D:\Github\TES-Programs\Elmer-Projects\results\case_phase24_short_bdf2_reuse_off_1us\case_phase24_short_bdf2_reuse_off_1us_series.csv`
- Reference: `CPU/MUMPS Phase23` — `D:\Github\TES-Programs\Elmer-Projects\results\case_p19_pulse_phase23_tight\case_p19_pulse_phase23_tight_series.csv`
- Event: `0.02002 s`
- Common post-event window: `0..0.901 µs` (20 rows)
- Baselines: candidate `203.747197482 µA`, reference `147.375245042 µA`, offset `+56.371952440 µA`

## Difference onset

- First raw absolute difference ≥ `0.1 µA`: **3.63797881e-12 µs**
- First baseline-corrected absolute difference ≥ `0.1 µA`: **3.63797881e-12 µs**
- First raw relative difference ≥ `0.1%`: **3.63797881e-12 µs**
- Maximum raw difference: `+58.874304002 µA` at `3.63797881e-12 µs`
- Maximum baseline-corrected difference: `+2.502351563 µA` at `3.63797881e-12 µs`

## Checkpoints

| t from pulse (µs) | candidate (µA) | reference (µA) | raw Δ (µA) | corrected Δ (µA) |
|---:|---:|---:|---:|---:|
| 0.1 | 206.215795451 | 147.389028372 | +58.826767079 | +2.454814640 |
| 0.2 | 206.099097185 | 147.386253795 | +58.712843390 | +2.340890951 |
| 0.3 | 205.936803193 | 147.380147953 | +58.556655241 | +2.184702801 |
| 0.4 | 205.747580727 | 147.370795763 | +58.376784964 | +2.004832525 |
| 0.5 | 205.523954393 | 147.358717780 | +58.165236614 | +1.793284174 |
| 0.6 | 205.289737268 | 147.344626368 | +57.945110900 | +1.573158461 |
| 0.7 | 205.020605108 | 147.329182491 | +57.691422617 | +1.319470177 |
| 0.8 | 204.746003045 | 147.312935484 | +57.433067561 | +1.061115122 |
| 0.9 | 204.459190896 | 147.296312315 | +57.162878581 | +0.790926141 |

The raw difference includes the pre-pulse operating-point offset. The corrected difference subtracts each run's own mean over the configured pre-pulse baseline window.
The comparison is based on recorded output samples; it does not treat internal nonlinear/Krylov epochs as physical time steps.

Files: `summary.json`, `aligned_difference.csv`, `checkpoints.csv`, `comparison.png`
