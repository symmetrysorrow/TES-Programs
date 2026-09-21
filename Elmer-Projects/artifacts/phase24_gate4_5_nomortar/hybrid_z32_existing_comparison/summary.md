# Optimized SinglePixel existing hybrid Stycast32 mortar HYPRE/MUMPS versus COMSOL

- Comparison window: 0–100.000 µs after the 20.020 ms pulse
- Mesh: `mesh_singlepixel_prod_v2` (optimized production-v2)
- Solver: `existing hybrid Stycast32 mortar HYPRE/MUMPS`
- Early timestep: 0.625 µs (optimized hybrid grid)
- COMSOL baseline: 143.055049 µA
- existing hybrid Stycast32 mortar HYPRE/MUMPS baseline: 147.355140 µA (+3.006%)
- Maximum absolute waveform difference: 0.035783 µA at 98.750 µs (0.460% of COMSOL full-trace peak)
- RMSE: 0.027995 µA (0.360% of COMSOL full-trace peak)
- t10 (COMSOL / existing hybrid Stycast32 mortar HYPRE/MUMPS): 41.7428 / 41.1820 µs
- t50 (COMSOL / existing hybrid Stycast32 mortar HYPRE/MUMPS): 92.8374 / 92.1898 µs

The traces are compared after subtracting each model's own pre-pulse baseline.
