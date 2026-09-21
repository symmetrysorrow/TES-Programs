# Optimized SinglePixel Phase24 HYPRE 1ms versus COMSOL

- Comparison window: 0–1000.000 µs after the 20.020 ms pulse
- Mesh: `mesh_singlepixel_prod_v2` (optimized production-v2)
- Solver: `Phase24 HYPRE 1ms`
- Early timestep: 0.625 µs (optimized hybrid grid)
- COMSOL baseline: 143.055049 µA
- Phase24 HYPRE 1ms baseline: 147.801439 µA (+3.318%)
- Maximum absolute waveform difference: 25.263894 µA at 1000.000 µs (324.944% of COMSOL full-trace peak)
- RMSE: 14.236803 µA (183.114% of COMSOL full-trace peak)
- t10 (COMSOL / Phase24 HYPRE 1ms): 41.7428 / 0.6441 µs
- t50 (COMSOL / Phase24 HYPRE 1ms): 92.8374 / 2.8282 µs

The traces are compared after subtracting each model's own pre-pulse baseline.
