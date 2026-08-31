# Optimized SinglePixel AMGX versus COMSOL

- Comparison window: 0–100.000 µs after the 20.020 ms pulse
- Mesh: `mesh_singlepixel_prod_v2` (optimized production-v2)
- Early timestep: 0.625 µs (optimized hybrid grid)
- COMSOL baseline: 143.055049 µA
- AMGX baseline: 143.687322 µA (+0.442%)
- Maximum absolute waveform difference: 3.984134 µA at 100.000 µs (51.244% of COMSOL full-trace peak)
- RMSE: 1.942205 µA (24.981% of COMSOL full-trace peak)
- t10 (COMSOL / AMGX): 41.7428 / n/a µs
- t50 (COMSOL / AMGX): 92.8374 / n/a µs

The traces are compared after subtracting each model's own pre-pulse baseline.
