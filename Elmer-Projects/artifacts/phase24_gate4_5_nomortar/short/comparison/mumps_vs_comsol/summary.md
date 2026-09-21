# Optimized SinglePixel no-mortar MUMPS versus COMSOL

- Comparison window: 0–0.900 µs after the 20.020 ms pulse
- Mesh: `mesh_singlepixel_prod_v2` (optimized production-v2)
- Solver: `no-mortar MUMPS`
- Early timestep: 0.625 µs (optimized hybrid grid)
- COMSOL baseline: 143.055049 µA
- no-mortar MUMPS baseline: 144.742692 µA (+1.180%)
- Maximum absolute waveform difference: 0.080698 µA at 0.900 µs (1.038% of COMSOL full-trace peak)
- RMSE: 0.038556 µA (0.496% of COMSOL full-trace peak)
- t10 (COMSOL / no-mortar MUMPS): 41.7428 / n/a µs
- t50 (COMSOL / no-mortar MUMPS): 92.8374 / n/a µs

The traces are compared after subtracting each model's own pre-pulse baseline.
