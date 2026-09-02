# Optimized SinglePixel GPU AMGX penalty1e8 vs CPU Phase19 versus COMSOL

- Comparison window: 0–100.000 µs after the 20.020 ms pulse
- Mesh: `mesh_singlepixel_prod_v2` (optimized production-v2)
- Solver: `GPU AMGX penalty1e8 vs CPU Phase19`
- Early timestep: 0.625 µs (optimized hybrid grid)
- COMSOL baseline: 147.375617 µA
- GPU AMGX penalty1e8 vs CPU Phase19 baseline: 147.357719 µA (-0.012%)
- Maximum absolute waveform difference: 4.450717 µA at 100.000 µs (56.471% of COMSOL full-trace peak)
- RMSE: 2.952052 µA (37.456% of COMSOL full-trace peak)
- t10 (COMSOL / GPU AMGX penalty1e8 vs CPU Phase19): 6.4732 / 0.2058 µs
- t50 (COMSOL / GPU AMGX penalty1e8 vs CPU Phase19): 42.1462 / n/a µs

The traces are compared after subtracting each model's own pre-pulse baseline.
