# Optimized SinglePixel nonconforming Stycast32 mortar MUMPS versus COMSOL

- Comparison window: 0–38.750 µs after the 20.020 ms pulse
- Mesh: `mesh_singlepixel_prod_v2` (optimized production-v2)
- Solver: `nonconforming Stycast32 mortar MUMPS`
- Early timestep: 0.625 µs (optimized hybrid grid)
- COMSOL baseline: 143.055049 µA
- nonconforming Stycast32 mortar MUMPS baseline: 165.706588 µA (+15.834%)
- Maximum absolute waveform difference: 17.826591 µA at 38.750 µs (229.286% of COMSOL full-trace peak)
- RMSE: 12.445186 µA (160.070% of COMSOL full-trace peak)
- t10 (COMSOL / nonconforming Stycast32 mortar MUMPS): 41.7428 / n/a µs
- t50 (COMSOL / nonconforming Stycast32 mortar MUMPS): 92.8374 / n/a µs

The traces are compared after subtracting each model's own pre-pulse baseline.
