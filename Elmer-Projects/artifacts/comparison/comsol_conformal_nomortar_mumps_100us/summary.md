# Single-pixel TES current comparison

Reference: COMSOL `docs/Single-Pixel.txt`; Elmer: `results\case_tes_pulse_singlepixel_conformal_gpu_fine_hybrid_177step\case_tes_pulse_singlepixel_conformal_gpu_fine_hybrid_177step_series.csv`.

| Metric | COMSOL | Elmer | Elmer error vs COMSOL |
|---|---:|---:|---:|
| Baseline current [µA] | 143.055049 | 154.851455 | +8.25% |
| Minimum current [µA] | 135.280206 | 148.585783 | +9.84% |
| Peak current drop [µA] | 7.774844 | 6.265672 | -19.41% |
| Peak delay from pulse [ms] | 0.428000 | 0.100001 | -76.64% |
| 10–90% rise time [ms] | 0.161440 | 0.074037 | -54.14% |

Baseline window: 19.50–20.02 ms. The peak is the maximum post-pulse current decrease; 10%/90% crossings are linearly interpolated.
