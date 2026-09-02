# Optimized SinglePixel AMGX versus COMSOL

- Comparison window: 0–100.000 µs after the 20.020 ms pulse
- Mesh: `mesh_singlepixel_prod_v2` (optimized production-v2)
- Early timestep: 0.625 µs (optimized hybrid grid)
- COMSOL baseline: 143.055049 µA
- AMGX baseline: 191.669672 µA (+33.983%)
- Maximum absolute waveform difference: 2.794781 µA at 100.000 µs (35.946% of COMSOL full-trace peak)
- RMSE: 1.190085 µA (15.307% of COMSOL full-trace peak)
- t10 (COMSOL / AMGX): 41.7428 / 37.2630 µs
- t50 (COMSOL / AMGX): 92.8374 / n/a µs

The traces are compared after subtracting each model's own pre-pulse baseline.
