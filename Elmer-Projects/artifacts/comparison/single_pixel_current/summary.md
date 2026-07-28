# Single-pixel TES current comparison

Reference: COMSOL `docs/Single-Pixel.txt`; Elmer: `artifacts\series\tes_pulse_20ms_3x_inner_circuit_pulse_aligned_partial.csv`.

| Metric | COMSOL | Elmer | Elmer error vs COMSOL |
|---|---:|---:|---:|
| Baseline current [µA] | 143.055049 | 148.150677 | +3.56% |
| Minimum current [µA] | 135.280206 | 148.149300 | +9.51% |
| Peak current drop [µA] | 7.774844 | 0.001377 | -99.98% |
| Peak delay from pulse [ms] | 0.428000 | 2.130000 | +397.66% |
| 10–90% rise time [ms] | 0.161440 | 2.116231 | +1210.84% |

Baseline window: 19.50–20.02 ms. The peak is the maximum post-pulse current decrease; 10%/90% crossings are linearly interpolated.
