# Optimized SinglePixel production-v2 direct MUMPS versus COMSOL

- Comparison window: 0–0.900 µs after the 20.020 ms pulse
- Mesh: `mesh_singlepixel_prod_v2` (optimized production-v2)
- Solver: `production-v2 direct MUMPS`
- Early timestep: 0.625 µs (optimized hybrid grid)
- COMSOL baseline: 143.055049 µA
- production-v2 direct MUMPS baseline: 143.777839 µA (+0.505%)
- Maximum absolute waveform difference: 0.000437 µA at 0.000 µs (0.006% of COMSOL full-trace peak)
- RMSE: 0.000329 µA (0.004% of COMSOL full-trace peak)
- t10 (COMSOL / production-v2 direct MUMPS): 41.7428 / n/a µs
- t50 (COMSOL / production-v2 direct MUMPS): 92.8374 / n/a µs

The traces are compared after subtracting each model's own pre-pulse baseline.
