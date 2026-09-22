# Historical mesh strict steady validation

- Historical saved restart: **168.563175809 mK / 143.774677487 uA**.
- Restart integrity gate: **PASS**; load, circuit-init, and pre-first-assembly fields agree.
- Initial captured full residual: L2=6.332096599e-10; primal=6.332096599e-10; constraint=6.202569190e-20; max=1.699747543e-10.
- One-shot thermal correction: Delta T=-0.000010701 mK; field-revaluated Delta I=0.000372293 uA; post-field current=143.775049779 uA.
- Full native MUMPS refinement: 168.563143389 mK / 143.775805366 uA after 21 nonlinear rows (interrupted_after_capture; no final result file).
- Final captured residual: L2=6.332090280e-10; primal=6.332090280e-10; constraint=3.328162390e-20; max=1.699745613e-10.

The historical 143.78-uA restart is a strict fixed-point candidate under the current native HeatSolve + CPU MUMPS replay: the one-shot correction is only 1.07e-5 mK and 3.72e-4 uA, and 21 captured nonlinear rows remain in the same 143.776-uA basin. The pre-solve residual is larger than the direct-solve residual because the saved restart is not the exact linear-system solution, but its physical correction is negligible.
- FE boundary diagnostic: bath out=3.205893940e-10 W historical vs 1.130754445e-10 W Phase24; secant G_eff=1.727019030e-08 vs 6.829454428e-09 W/K.
- One-sided historical_refined_or_one_shot TES_to_membrane flux: 8.602378947e-10 W / -3.596213918e-10 W; signed mismatch=5.006165029e-10 W (diagnostic wedge/mortar extraction, not solver-native flux).
- One-sided historical_refined_or_one_shot TES_to_Stycast flux: 1.391377180e-09 W / 8.782559450e-16 W; signed mismatch=1.391378058e-09 W (diagnostic wedge/mortar extraction, not solver-native flux).
- One-sided phase24_refined TES_to_membrane flux: 1.286061582e-09 W / -2.317331670e-10 W; signed mismatch=1.054328415e-09 W (diagnostic wedge/mortar extraction, not solver-native flux).
- One-sided phase24_refined TES_to_Stycast flux: 1.113879811e-09 W / 6.211380242e-16 W; signed mismatch=1.113880432e-09 W (diagnostic wedge/mortar extraction, not solver-native flux).

HYPRE/GPU tuning remains NO-GO.
