# Phase24 difference evidence

Generated (UTC): `2026-09-15T17:53:24.096410+00:00`
Git revision: `45778dd1f918b861eb6b7b342672336f7c28b5bd`

- Candidate: `same Phase24 / direct MUMPS` — `D:\Github\TES-Programs\Elmer-Projects\results\case_phase24_short_same_path_mumps_1us\case_phase24_short_same_path_mumps_1us_series.csv`
- Reference: `reference` — `D:\Github\TES-Programs\Elmer-Projects\results\case_p19_pulse_phase23_tight\case_p19_pulse_phase23_tight_series.csv`
- Event: `0.02002 s`
- Common post-event window: `0..0.901 µs` (20 rows)
- Baselines: candidate `144.268506298 µA`, reference `147.375245042 µA`, offset `-3.106738744 µA`

## Difference onset

- First raw absolute difference ≥ `0.1 µA`: **3.63797881e-12 µs**
- First baseline-corrected absolute difference ≥ `0.1 µA`: **3.63797881e-12 µs**
- First raw relative difference ≥ `0.1%`: **3.63797881e-12 µs**
- Maximum raw difference: `-4.072747815 µA` at `0.901 µs`
- Maximum baseline-corrected difference: `-0.966009071 µA` at `0.901 µs`

## Checkpoints

| t from pulse (µs) | candidate (µA) | reference (µA) | raw Δ (µA) | corrected Δ (µA) |
|---:|---:|---:|---:|---:|
| 0.1 | 144.053174362 | 147.389028372 | -3.335854010 | -0.229115266 |
| 0.2 | 144.029251023 | 147.386253795 | -3.357002772 | -0.250264028 |
| 0.3 | 143.981552123 | 147.380147953 | -3.398595829 | -0.291857086 |
| 0.4 | 143.904187162 | 147.370795763 | -3.466608601 | -0.359869857 |
| 0.5 | 143.798675385 | 147.358717780 | -3.560042394 | -0.453303650 |
| 0.6 | 143.671065493 | 147.344626368 | -3.673560875 | -0.566822131 |
| 0.7 | 143.528744995 | 147.329182491 | -3.800437496 | -0.693698752 |
| 0.8 | 143.378354761 | 147.312935484 | -3.934580723 | -0.827841979 |
| 0.9 | 143.224932750 | 147.296312315 | -4.071379566 | -0.964640822 |

The raw difference includes the pre-pulse operating-point offset. The corrected difference subtracts each run's own mean over the configured pre-pulse baseline window.
The comparison is based on recorded output samples; it does not treat internal nonlinear/Krylov epochs as physical time steps.

Files: `summary.json`, `aligned_difference.csv`, `checkpoints.csv`, `comparison.png`
