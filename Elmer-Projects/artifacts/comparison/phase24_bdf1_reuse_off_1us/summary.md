# Phase24 difference evidence

Generated (UTC): `2026-09-15T16:30:25.066007+00:00`
Git revision: `a2bfedc466662a8a68598f95e960bf6bcf224942`

- Candidate: `BDF1 / reuse OFF` — `D:\Github\TES-Programs\Elmer-Projects\results\case_phase24_short_bdf1_reuse_off_1us\case_phase24_short_bdf1_reuse_off_1us_series.csv`
- Reference: `CPU/MUMPS Phase23` — `D:\Github\TES-Programs\Elmer-Projects\results\case_p19_pulse_phase23_tight\case_p19_pulse_phase23_tight_series.csv`
- Event: `0.02002 s`
- Common post-event window: `0..0.901 µs` (20 rows)
- Baselines: candidate `203.932518710 µA`, reference `147.375245042 µA`, offset `+56.557273668 µA`

## Difference onset

- First raw absolute difference ≥ `0.1 µA`: **3.63797881e-12 µs**
- First baseline-corrected absolute difference ≥ `0.1 µA`: **3.63797881e-12 µs**
- First raw relative difference ≥ `0.1%`: **3.63797881e-12 µs**
- Maximum raw difference: `+59.756572515 µA` at `0.021 µs`
- Maximum baseline-corrected difference: `+3.199298847 µA` at `0.021 µs`

## Checkpoints

| t from pulse (µs) | candidate (µA) | reference (µA) | raw Δ (µA) | corrected Δ (µA) |
|---:|---:|---:|---:|---:|
| 0.1 | 207.059015382 | 147.389028372 | +59.669987010 | +3.112713342 |
| 0.2 | 206.926937894 | 147.386253795 | +59.540684100 | +2.983410431 |
| 0.3 | 206.798613335 | 147.380147953 | +59.418465382 | +2.861191714 |
| 0.4 | 206.676068381 | 147.370795763 | +59.305272618 | +2.747998950 |
| 0.5 | 206.546643736 | 147.358717780 | +59.187925956 | +2.630652288 |
| 0.6 | 206.403730755 | 147.344626368 | +59.059104387 | +2.501830719 |
| 0.7 | 206.279746246 | 147.329182491 | +58.950563755 | +2.393290086 |
| 0.8 | 206.128401184 | 147.312935484 | +58.815465701 | +2.258192032 |
| 0.9 | 205.999050425 | 147.296312315 | +58.702738110 | +2.145464442 |

The raw difference includes the pre-pulse operating-point offset. The corrected difference subtracts each run's own mean over the configured pre-pulse baseline window.
The comparison is based on recorded output samples; it does not treat internal nonlinear/Krylov epochs as physical time steps.

Files: `summary.json`, `aligned_difference.csv`, `checkpoints.csv`, `comparison.png`
