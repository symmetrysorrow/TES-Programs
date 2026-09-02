# Optimized SinglePixel AMGX versus COMSOL

- Comparison window: 0–100.000 µs after the 20.020 ms pulse
- Mesh: `mesh_singlepixel_conformal_gpu_fine` (90,872 nodes, all-tetra, shared contact nodes)
- Time grid: optimized 177-step hybrid grid (0.625 µs fine stage)
- COMSOL baseline: 143.055049 µA
- AMGX baseline: 154.851462 µA (+8.246%)
- Maximum absolute waveform difference: 4.126886 µA at 100.000 µs (53.080% of COMSOL full-trace peak)
- RMSE: 2.037903 µA (26.212% of COMSOL full-trace peak)
- t10 (COMSOL / AMGX): 41.7428 / n/a µs
- t50 (COMSOL / AMGX): 92.8374 / n/a µs

The traces are compared after subtracting each model's own pre-pulse baseline.
