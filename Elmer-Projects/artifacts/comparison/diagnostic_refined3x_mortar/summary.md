# Single-pixel TES current comparison

Reference: COMSOL `docs/Single-Pixel.txt`; Elmer: `results\case_tes_pulse_20ms_3x_comsol_grid\tes_pulse_20ms_3x_comsol_grid_series.csv`.

| Metric | COMSOL | Elmer | Elmer error vs COMSOL |
|---|---:|---:|---:|
| Baseline current [µA] | 143.055049 | 148.159981 | +3.57% |
| Minimum current [µA] | 135.280206 | 148.159616 | +9.52% |
| Peak current drop [µA] | 7.774844 | 0.000364 | -100.00% |
| Peak delay from pulse [ms] | 0.428000 | 0.000012 | -100.00% |
| 10–90% rise time [ms] | 0.161440 | 0.002147 | -98.67% |

Baseline window: 19.50–20.02 ms. The peak is the maximum post-pulse current decrease; 10%/90% crossings are linearly interpolated.
