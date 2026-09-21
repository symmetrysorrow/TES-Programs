# Optimized SinglePixel Phase24 HYPRE forced AMG/Krylov setup versus COMSOL

- Comparison window: 0–0.900 µs after the 20.020 ms pulse
- Mesh: `mesh_singlepixel_prod_v2` (optimized production-v2)
- Solver: `Phase24 HYPRE forced AMG/Krylov setup`
- Early timestep: 0.625 µs (optimized hybrid grid)
- COMSOL baseline: 143.055049 µA
- Phase24 HYPRE forced AMG/Krylov setup baseline: 147.787052 µA (+3.308%)
- Maximum absolute waveform difference: 1.223343 µA at 0.900 µs (15.735% of COMSOL full-trace peak)
- RMSE: 0.690717 µA (8.884% of COMSOL full-trace peak)
- t10 (COMSOL / Phase24 HYPRE forced AMG/Krylov setup): 41.7428 / 0.6189 µs
- t50 (COMSOL / Phase24 HYPRE forced AMG/Krylov setup): 92.8374 / n/a µs

The traces are compared after subtracting each model's own pre-pulse baseline.
