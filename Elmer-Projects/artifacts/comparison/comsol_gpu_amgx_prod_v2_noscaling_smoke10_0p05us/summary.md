# Optimized SinglePixel AMGX / RTX 3060 Ti / no-scaling smoke10 versus COMSOL

- Comparison window: 0–0.050 µs after the 20.020 ms pulse
- Mesh: `mesh_singlepixel_prod_v2` (optimized production-v2)
- Solver: `AMGX / RTX 3060 Ti / no-scaling smoke10`
- Early timestep: 0.625 µs (optimized hybrid grid)
- COMSOL baseline: 143.055049 µA
- AMGX / RTX 3060 Ti / no-scaling smoke10 baseline: 143.774169 µA (+0.503%)
- Maximum absolute waveform difference: 0.264739 µA at 0.050 µs (3.405% of COMSOL full-trace peak)
- RMSE: 0.187200 µA (2.408% of COMSOL full-trace peak)
- t10 (COMSOL / AMGX / RTX 3060 Ti / no-scaling smoke10): 41.7428 / n/a µs
- t50 (COMSOL / AMGX / RTX 3060 Ti / no-scaling smoke10): 92.8374 / n/a µs

The traces are compared after subtracting each model's own pre-pulse baseline.
