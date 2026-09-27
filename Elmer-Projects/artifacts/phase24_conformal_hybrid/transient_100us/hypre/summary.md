# Optimized SinglePixel hypre versus COMSOL

- Comparison window: 0–99.376 µs after the 20.020 ms pulse
- Mesh: `mesh_singlepixel_prod_v2` (optimized production-v2)
- Solver: `hypre`
- Early timestep: 0.625 µs (optimized hybrid grid)
- COMSOL baseline: 143.055049 µA
- hypre baseline: 139.676818 µA (-2.361%)
- Maximum absolute waveform difference: 0.080263 µA at 99.376 µs (1.032% of COMSOL full-trace peak)
- RMSE: 0.025834 µA (0.332% of COMSOL full-trace peak)
- t10 (COMSOL / hypre): 41.7428 / 41.5138 µs
- t50 (COMSOL / hypre): 92.8374 / 93.8182 µs

The traces are compared after subtracting each model's own pre-pulse baseline.
