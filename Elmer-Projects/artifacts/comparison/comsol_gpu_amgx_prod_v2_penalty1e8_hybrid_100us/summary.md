# Optimized SinglePixel AMGX / RTX 3060 Ti / production-v2 penalty1e8 versus COMSOL

- Comparison window: 0–100.000 µs after the 20.020 ms pulse
- Mesh: `mesh_singlepixel_prod_v2` (optimized production-v2)
- Solver: `AMGX / RTX 3060 Ti / production-v2 penalty1e8`
- Early timestep: 0.625 µs (optimized hybrid grid)
- COMSOL baseline: 143.055049 µA
- AMGX / RTX 3060 Ti / production-v2 penalty1e8 baseline: 143.777029 µA (+0.505%)
- Maximum absolute waveform difference: 5.313753 µA at 100.000 µs (68.345% of COMSOL full-trace peak)
- RMSE: 2.988744 µA (38.441% of COMSOL full-trace peak)
- t10 (COMSOL / AMGX / RTX 3060 Ti / production-v2 penalty1e8): 41.7428 / n/a µs
- t50 (COMSOL / AMGX / RTX 3060 Ti / production-v2 penalty1e8): 92.8374 / n/a µs

The traces are compared after subtracting each model's own pre-pulse baseline.
