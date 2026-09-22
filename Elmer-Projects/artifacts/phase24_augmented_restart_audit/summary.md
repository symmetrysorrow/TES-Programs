# Phase24 augmented restart audit

The corrected Gate3 steady restart loads the circuit checkpoint successfully before the first HeatSolve assembly. The first mismatch is none for the persisted electrical fields: Current, PreviousCurrent, Resistance, and Power are identical at the saved checkpoint, Stage A, Stage B, and Stage C.

| field | saved | load_after | circuit_init_after | pre_first_assembly |
|---|---:|---:|---:|---:|
| Current (uA) | 143.568114993 | 143.568114993 | 143.568114993 | 143.568114993 |
| PreviousCurrent (uA) | 143.568114993 | 143.568114993 | 143.568114993 | 143.568114993 |
| Resistance (ohm) | 0.015522836332 | 0.015522836332 | 0.015522836332 | 0.015522836332 |
| Power (W) | 3.199063314934e-10 | 3.199063314934e-10 | 3.199063314934e-10 | 3.199063314934e-10 |
| TES T (mK) | 168.569114993 | 168.569114993 | 168.569114993 | 168.569114993 |

Runtime reported `file_exists=T`, `open_success=T`, `parse_success=T`, and `loaded=T`. The configured path is repository-relative and resolves from ElmerSolver cwd to the isolated diagnostic snapshot. The original Gate3 state SHA256 remained `9c2eca9bbb26fc6fab1949511787b87350ccdc370658ed2a41d458959f5971a6`; it was not overwritten.

Root cause had two parts: the diagnostic harness rewrote `work/meshes/...` to `../../work/...`, and native `HeatSolve.F90::TESInnerCircuitUpdate` only read state when `TransientSimulation` was true. The minimal fixes preserve the repository-relative path and allow an explicitly supplied checkpoint to be read by a steady restart; fresh steady cases without a state file still use the existing default initialization.

The corrected Gate3 one-shot direct image is `Delta T = -6.199616730 mK` and field-re-evaluated `Delta I = +246.593147 uA`, with full captured solve residual L2 `6.86e-16`. The old one-shot numbers are superseded.

The corrected refinement began from the saved 218.667665604 uA state, produced 12 captured circuit rows, and ended at 166.557530019 mK / 218.643409371 uA in the last captured iterate. The native captured full residual was L2 `6.19e-16`, max `1.15e-16`; constraint residual L2 was `2.45e-22`. This is a single corrected basin; it does not establish a second steady branch.

No physics parameters, materials, TES law, mesh, mortar formulation, HYPRE/GPU settings, or production timestep were changed. HYPRE/GPU tuning should remain paused until this corrected refinement is accepted as the restart baseline. The final-series issue now has a launcher-side normal-completion-only flush: it appends the final iteration row only when newer than the existing series tail, and leaves failed runs and already-flushed series unchanged.
