# Phase24 difference evidence

Generated (UTC): `2026-09-15T17:53:25.126010+00:00`
Git revision: `45778dd1f918b861eb6b7b342672336f7c28b5bd`

- Candidate: `same Phase24 / direct MUMPS` — `D:\Github\TES-Programs\Elmer-Projects\results\case_phase24_short_same_path_mumps_1us\case_phase24_short_same_path_mumps_1us_series.csv`
- Reference: `reference` — `D:\Github\TES-Programs\Elmer-Projects\results\case_phase24_short_bdf2_reuse_on_1us\case_phase24_short_bdf2_reuse_on_1us_series.csv`
- Event: `0.02002 s`
- Common post-event window: `0..0.901 µs` (20 rows)
- Baselines: candidate `144.268506298 µA`, reference `147.801439420 µA`, offset `-3.532933122 µA`

## Difference onset

- First raw absolute difference ≥ `0.1 µA`: **3.63797881e-12 µs**
- First baseline-corrected absolute difference ≥ `0.1 µA`: **3.63797881e-12 µs**
- First raw relative difference ≥ `0.1%`: **3.63797881e-12 µs**
- Maximum raw difference: `-3.446224960 µA` at `0.201 µs`
- Maximum baseline-corrected difference: `+0.142527967 µA` at `0.901 µs`

## Checkpoints

| t from pulse (µs) | candidate (µA) | reference (µA) | raw Δ (µA) | corrected Δ (µA) |
|---:|---:|---:|---:|---:|
| 0.1 | 144.053174362 | 147.496189303 | -3.443014941 | +0.089918180 |
| 0.2 | 144.029251023 | 147.475445111 | -3.446194088 | +0.086739034 |
| 0.3 | 143.981552123 | 147.425217760 | -3.443665637 | +0.089267485 |
| 0.4 | 143.904187162 | 147.341609721 | -3.437422559 | +0.095510563 |
| 0.5 | 143.798675385 | 147.227404131 | -3.428728746 | +0.104204376 |
| 0.6 | 143.671065493 | 147.090694970 | -3.419629477 | +0.113303645 |
| 0.7 | 143.528744995 | 146.939239590 | -3.410494595 | +0.122438527 |
| 0.8 | 143.378354761 | 146.781606123 | -3.403251362 | +0.129681760 |
| 0.9 | 143.224932750 | 146.615465644 | -3.390532894 | +0.142400227 |

The raw difference includes the pre-pulse operating-point offset. The corrected difference subtracts each run's own mean over the configured pre-pulse baseline window.
The comparison is based on recorded output samples; it does not treat internal nonlinear/Krylov epochs as physical time steps.

Files: `summary.json`, `aligned_difference.csv`, `checkpoints.csv`, `comparison.png`
