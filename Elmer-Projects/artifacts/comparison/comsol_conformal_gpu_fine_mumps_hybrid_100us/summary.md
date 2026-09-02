# Optimized SinglePixel GPU binary / MUMPS versus COMSOL

- Comparison window: 0–100.000 µs after the 20.020 ms pulse
- Mesh: `mesh_singlepixel_prod_v2` (optimized production-v2)
- Solver: `GPU binary / MUMPS`
- Early timestep: 0.625 µs (optimized hybrid grid)
- COMSOL baseline: 143.055049 µA
- GPU binary / MUMPS baseline: 154.851455 µA (+8.246%)
- Maximum absolute waveform difference: 2.938712 µA at 48.000 µs (37.798% of COMSOL full-trace peak)
- RMSE: 2.372670 µA (30.517% of COMSOL full-trace peak)
- t10 (COMSOL / GPU binary / MUMPS): 41.7428 / 7.8348 µs
- t50 (COMSOL / GPU binary / MUMPS): 92.8374 / 44.9400 µs

The traces are compared after subtracting each model's own pre-pulse baseline.
