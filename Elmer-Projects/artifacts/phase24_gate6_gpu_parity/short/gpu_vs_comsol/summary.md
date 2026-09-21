# Optimized SinglePixel HIP HYPRE GPU versus COMSOL

- Comparison window: 0–0.900 µs after the 20.020 ms pulse
- Mesh: `mesh_singlepixel_prod_v2` (optimized production-v2)
- Solver: `HIP HYPRE GPU`
- Early timestep: 0.625 µs (optimized hybrid grid)
- COMSOL baseline: 143.055049 µA
- HIP HYPRE GPU baseline: 144.742689 µA (+1.180%)
- Maximum absolute waveform difference: 0.001616 µA at 0.900 µs (0.021% of COMSOL full-trace peak)
- RMSE: 0.001050 µA (0.014% of COMSOL full-trace peak)
- t10 (COMSOL / HIP HYPRE GPU): 41.7428 / n/a µs
- t50 (COMSOL / HIP HYPRE GPU): 92.8374 / n/a µs

The traces are compared after subtracting each model's own pre-pulse baseline.
