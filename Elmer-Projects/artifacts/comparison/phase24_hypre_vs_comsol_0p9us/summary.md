# Optimized SinglePixel Phase24 HYPRE versus COMSOL

- Comparison window: 0–0.900 µs after the 20.020 ms pulse
- Mesh: `mesh_singlepixel_prod_v2` (optimized production-v2)
- Solver: `Phase24 HYPRE`
- Early timestep: 0.625 µs (optimized hybrid grid)
- COMSOL baseline: 143.055049 µA
- Phase24 HYPRE baseline: 147.801439 µA (+3.318%)
- Maximum absolute waveform difference: 1.185975 µA at 0.900 µs (15.254% of COMSOL full-trace peak)
- RMSE: 0.671649 µA (8.639% of COMSOL full-trace peak)
- t10 (COMSOL / Phase24 HYPRE): 41.7428 / 0.6441 µs
- t50 (COMSOL / Phase24 HYPRE): 92.8374 / n/a µs

The traces are compared after subtracting each model's own pre-pulse baseline.
