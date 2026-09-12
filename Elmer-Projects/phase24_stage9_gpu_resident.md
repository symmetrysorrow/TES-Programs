# Phase24 Stage 9: persistent GPU HYPRE path

Stage 9 keeps FEM assembly on the CPU and moves the dominant linear-algebra work into a persistent HYPRE GPU execution path. CUDA and HIP are selected at build time from HYPRE's backend macros; a GPU request compiled against a CPU-only HYPRE library is converted to the existing CPU path instead of failing.

The Phase24 container now retains the IJ matrix, ParCSR matrix, RHS and warm-start solution vectors, FlexGMRES state, BoomerAMG state, and device work allocations across solves. Matrix and RHS changes are inserted through the host boundary, then HYPRE executes the solve and preconditioner work with the selected device. The normal path does not recreate the HYPRE objects. Case A, B1, B2, recovery, epoch commits, lagged AMG, and warm starts remain in the same lifecycle as the optimized CPU path.

The backend runner is `scripts/run_phase24_hypre_wsl.ps1`. It accepts `cpu`, `cuda`, or `hip`, keeps HYPRE installations side-by-side by release, and defaults to the validated CUDA HYPRE 3.0.0 build. The generated GPU case sets `HYPRE GPU = True` and uses FlexGMRES with BoomerAMG.

## Validation

The seven-step CUDA production run completed on the NVIDIA GeForce RTX 3060 Ti with exit code 0 and `ALL DONE`. It produced the same final result norm as the optimized CPU run to a relative difference of `1.54e-12`, while wall time fell from 355.41 s to 193.62 s (1.84x speedup).

The complete counters, device confirmation, lifecycle counts, transfer-boundary counts, and CPU comparison are recorded in `phase24_stage9.json`.

HYPRE 3.1.0 was also tested. Its CUDA build reported the real device but aborted during the first sparse setup with `cudaErrorInvalidDevice`; production therefore uses the known-good 3.0.0 CUDA build by default, while 3.1.0 remains an explicit compatibility-test option. Explicit IJ migration is disabled for this reason; the validated path uses HYPRE's supported host-IJ insertion/device-execution boundary and retains the persistent device-side solver state.
