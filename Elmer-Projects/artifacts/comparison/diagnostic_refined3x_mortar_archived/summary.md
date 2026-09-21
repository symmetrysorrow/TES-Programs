# Single-pixel TES current comparison

Reference: COMSOL `docs/Single-Pixel.txt`; Elmer: `artifacts\series\tes_pulse_20ms_3x_series.csv`.

| Metric | COMSOL | Elmer | Elmer error vs COMSOL |
|---|---:|---:|---:|
| Baseline current [µA] | 143.055049 | 148.130158 | +3.55% |
| Minimum current [µA] | 135.280206 | 140.404956 | +3.79% |
| Peak current drop [µA] | 7.774844 | 7.725202 | -0.64% |
| Peak delay from pulse [ms] | 0.428000 | 0.500001 | +16.82% |
| 10–90% rise time [ms] | 0.161440 | 0.174026 | +7.80% |

Baseline window: 19.50–20.02 ms. The peak is the maximum post-pulse current decrease; 10%/90% crossings are linearly interpolated.
