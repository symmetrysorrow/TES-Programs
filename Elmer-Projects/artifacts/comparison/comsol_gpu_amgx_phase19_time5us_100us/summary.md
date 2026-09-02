# Optimized SinglePixel AMGX / RTX 3060 Ti / Phase19 hybrid-prism versus COMSOL

- Comparison window: 0–100.000 µs after the 20.020 ms pulse
- Mesh: `mesh_singlepixel_prod_v2` (optimized production-v2)
- Solver: `AMGX / RTX 3060 Ti / Phase19 hybrid-prism`
- Early timestep: 0.625 µs (optimized hybrid grid)
- COMSOL baseline: 143.055049 µA
- AMGX / RTX 3060 Ti / Phase19 hybrid-prism baseline: 147.383152 µA (+3.025%)
- Maximum absolute waveform difference: 11.202316 µA at 7.000 µs (144.084% of COMSOL full-trace peak)
- RMSE: 9.655376 µA (124.187% of COMSOL full-trace peak)
- t10 (COMSOL / AMGX / RTX 3060 Ti / Phase19 hybrid-prism): 41.7428 / 0.1941 µs
- t50 (COMSOL / AMGX / RTX 3060 Ti / Phase19 hybrid-prism): 92.8374 / 0.4193 µs

The traces are compared after subtracting each model's own pre-pulse baseline.
