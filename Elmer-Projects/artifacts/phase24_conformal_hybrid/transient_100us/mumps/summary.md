# Optimized SinglePixel mumps versus COMSOL

- Comparison window: 0–99.376 µs after the 20.020 ms pulse
- Mesh: `mesh_singlepixel_prod_v2` (optimized production-v2)
- Solver: `mumps`
- Early timestep: 0.625 µs (optimized hybrid grid)
- COMSOL baseline: 143.055049 µA
- mumps baseline: 139.676799 µA (-2.362%)
- Maximum absolute waveform difference: 0.080252 µA at 99.376 µs (1.032% of COMSOL full-trace peak)
- RMSE: 0.025829 µA (0.332% of COMSOL full-trace peak)
- t10 (COMSOL / mumps): 41.7428 / 41.5136 µs
- t50 (COMSOL / mumps): 92.8374 / 93.8180 µs

The traces are compared after subtracting each model's own pre-pulse baseline.
