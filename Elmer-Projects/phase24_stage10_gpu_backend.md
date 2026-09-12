# Phase24 Stage 10: HIP/ROCm production validation

Stage 10 brings the Radeon RX 9070 XT into the same persistent HYPRE lifecycle used by CPU and CUDA. FEM assembly remains CPU-side. Backend selection is below HeatSolve, so StructureEpoch, MatrixEpoch, HypreMatrixEpoch, PreconditionerEpoch, Case A/B1/B2/C, warm starts, adaptive AMG lagging, and recovery are shared.

## HIP validation

The installed AMD runtime is usable in WSL: `rocminfo` reports an AMD GPU agent named `AMD Radeon RX 9070 XT` with architecture `gfx1201`; HIP device allocation succeeds; and the rocSPARSE query passes. HYPRE 3.1.0 was built in an isolated prefix with HIP, rocSPARSE, rocRAND, rocBLAS, rocSOLVER, and rocTHRUST enabled.

The seven-step production run reports:

- backend: HIP
- device: AMD Radeon RX 9070 XT
- exit code: 0
- `ALL DONE`
- persistent GPU backend counter: 1
- matrix refreshes: 12
- FlexGMRES setups: 16
- real AMG setups: 10
- solves: 16
- recoveries: 3
- device allocations: 3

## Comparison

All backends produce the same final result to numerical precision. The CPU/GPU final norm relative difference is `1.54e-12`, with final relative change zero and the shared linear tolerance `5e-7`.

The conservative cold-start measurements are CUDA 193.62 s and HIP 197.93 s. After ROCm kernel caching, repeated production measurements are CUDA 169.83 s and HIP 169.51 s; CPU is 355.41 s. HIP is therefore selected as `HIP_BEST` for sustained repeated production, with CUDA retained as the explicit frozen reference and the better cold-start choice.

The complete measurements and evidence are in `phase24_stage10_gpu_backends.json`, `phase24_stage10_hip_validation.json`, and `phase24_stage10_backend_recommendation.json`.
