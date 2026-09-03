# Optimized SinglePixel CPU MUMPS versus COMSOL

- Comparison window: 0–179980.000 µs after the 20.020 ms pulse
- Mesh: `mesh_singlepixel_prod_v2` (optimized production-v2)
- Solver: `CPU MUMPS`
- Early timestep: 0.625 µs (optimized hybrid grid)
- COMSOL baseline: 143.055049 µA
- CPU MUMPS baseline: 143.777852 µA (+0.505%)
- Maximum absolute waveform difference: 0.750614 µA at 220.000 µs (9.654% of COMSOL full-trace peak)
- RMSE: 0.032824 µA (0.422% of COMSOL full-trace peak)
- t10 (COMSOL / CPU MUMPS): 41.7428 / 41.2819 µs
- t50 (COMSOL / CPU MUMPS): 92.8374 / 92.8493 µs

The traces are compared after subtracting each model's own pre-pulse baseline.
