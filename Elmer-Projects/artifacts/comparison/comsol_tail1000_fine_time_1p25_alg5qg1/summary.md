# Optimized SinglePixel CPU MUMPS (1.25 us prefix + validated tail) versus COMSOL

- Comparison window: 0–179980.000 µs after the 20.020 ms pulse
- Mesh: `mesh_singlepixel_prod_v2` (optimized production-v2)
- Solver: `CPU MUMPS (1.25 us prefix + validated tail)`
- Time grid: 0.625 us prefix + 1.25 us through 998.751 us + validated tail
- COMSOL baseline: 143.055049 µA
- CPU MUMPS (1.25 us prefix + validated tail) baseline: 143.774271 µA (+0.503%)
- Maximum absolute waveform difference: 0.065558 µA at 13980.000 µs (0.843% of COMSOL full-trace peak)
- RMSE: 0.021488 µA (0.276% of COMSOL full-trace peak)
- t10 (COMSOL / CPU MUMPS (1.25 us prefix + validated tail)): 41.7428 / 41.3509 µs
- t50 (COMSOL / CPU MUMPS (1.25 us prefix + validated tail)): 92.8374 / 92.9460 µs

The traces are compared after subtracting each model's own pre-pulse baseline.
