# Optimized SinglePixel Phase24 MUMPS versus COMSOL

- Comparison window: 0–0.900 µs after the 20.020 ms pulse
- Mesh: `mesh_singlepixel_prod_v2` (optimized production-v2)
- Solver: `Phase24 MUMPS`
- Early timestep: 0.625 µs (optimized hybrid grid)
- COMSOL baseline: 143.055049 µA
- Phase24 MUMPS baseline: 144.268506 µA (+0.848%)
- Maximum absolute waveform difference: 1.043575 µA at 0.900 µs (13.422% of COMSOL full-trace peak)
- RMSE: 0.568496 µA (7.312% of COMSOL full-trace peak)
- t10 (COMSOL / Phase24 MUMPS): 41.7428 / 0.7251 µs
- t50 (COMSOL / Phase24 MUMPS): 92.8374 / n/a µs

The traces are compared after subtracting each model's own pre-pulse baseline.
