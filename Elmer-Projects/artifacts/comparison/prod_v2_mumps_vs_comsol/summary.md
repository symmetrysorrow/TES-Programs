# Optimized SinglePixel Production v2 direct MUMPS versus COMSOL

- Comparison window: 0–100.000 µs after the 20.020 ms pulse
- Mesh: `mesh_singlepixel_prod_v2` (optimized production-v2)
- Solver: `Production v2 direct MUMPS`
- Early timestep: 0.625 µs (optimized hybrid grid)
- COMSOL baseline: 143.055049 µA
- Production v2 direct MUMPS baseline: 143.777756 µA (+0.505%)
- Maximum absolute waveform difference: 0.028629 µA at 36.000 µs (0.368% of COMSOL full-trace peak)
- RMSE: 0.015804 µA (0.203% of COMSOL full-trace peak)
- t10 (COMSOL / Production v2 direct MUMPS): 41.7428 / 41.2834 µs
- t50 (COMSOL / Production v2 direct MUMPS): 92.8374 / 92.8510 µs

The traces are compared after subtracting each model's own pre-pulse baseline.
