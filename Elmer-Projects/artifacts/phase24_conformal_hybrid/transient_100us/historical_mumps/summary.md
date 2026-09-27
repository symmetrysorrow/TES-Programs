# Optimized SinglePixel historical_mumps versus COMSOL

- Comparison window: 0–99.376 µs after the 20.020 ms pulse
- Mesh: `mesh_singlepixel_prod_v2` (optimized production-v2)
- Solver: `historical_mumps`
- Early timestep: 0.625 µs (optimized hybrid grid)
- COMSOL baseline: 143.055049 µA
- historical_mumps baseline: 143.777852 µA (+0.505%)
- Maximum absolute waveform difference: 0.028705 µA at 35.991 µs (0.369% of COMSOL full-trace peak)
- RMSE: 0.015913 µA (0.205% of COMSOL full-trace peak)
- t10 (COMSOL / historical_mumps): 41.7428 / 41.2819 µs
- t50 (COMSOL / historical_mumps): 92.8374 / 92.8493 µs

The traces are compared after subtracting each model's own pre-pulse baseline.
