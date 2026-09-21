# Optimized SinglePixel conformal Stycast32 stripped-internal-boundary HYPRE partial versus COMSOL

- Comparison window: 0–38.750 µs after the 20.020 ms pulse
- Mesh: `mesh_singlepixel_prod_v2` (optimized production-v2)
- Solver: `conformal Stycast32 stripped-internal-boundary HYPRE partial`
- Early timestep: 0.625 µs (optimized hybrid grid)
- COMSOL baseline: 143.055049 µA
- conformal Stycast32 stripped-internal-boundary HYPRE partial baseline: 143.534426 µA (+0.335%)
- Maximum absolute waveform difference: 0.617301 µA at 38.750 µs (7.940% of COMSOL full-trace peak)
- RMSE: 0.217332 µA (2.795% of COMSOL full-trace peak)
- t10 (COMSOL / conformal Stycast32 stripped-internal-boundary HYPRE partial): 41.7428 / n/a µs
- t50 (COMSOL / conformal Stycast32 stripped-internal-boundary HYPRE partial): 92.8374 / n/a µs

The traces are compared after subtracting each model's own pre-pulse baseline.
