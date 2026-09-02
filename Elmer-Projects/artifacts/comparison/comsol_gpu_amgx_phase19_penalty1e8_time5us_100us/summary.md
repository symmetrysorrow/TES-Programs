# Optimized SinglePixel AMGX / RTX 3060 Ti / Phase19 penalty1e8 versus COMSOL

- Comparison window: 0–100.000 µs after the 20.020 ms pulse
- Mesh: `mesh_singlepixel_prod_v2` (optimized production-v2)
- Solver: `AMGX / RTX 3060 Ti / Phase19 penalty1e8`
- Early timestep: 0.625 µs (optimized hybrid grid)
- COMSOL baseline: 143.055049 µA
- AMGX / RTX 3060 Ti / Phase19 penalty1e8 baseline: 147.357719 µA (+3.008%)
- Maximum absolute waveform difference: 2.283136 µA at 100.000 µs (29.366% of COMSOL full-trace peak)
- RMSE: 1.159239 µA (14.910% of COMSOL full-trace peak)
- t10 (COMSOL / AMGX / RTX 3060 Ti / Phase19 penalty1e8): 41.7428 / 0.2039 µs
- t50 (COMSOL / AMGX / RTX 3060 Ti / Phase19 penalty1e8): 92.8374 / n/a µs

The traces are compared after subtracting each model's own pre-pulse baseline.
