# Phase24 historical steady-parity audit

## Finding

The current full nonlinear Phase24 system has a residual-qualified solution at **218.643409 uA**.  The corrected Gate3 143.568115-uA checkpoint is not a fixed point: holding its accepted Joule source and solving the captured full constrained system changes the TES average by **-6.199616730 mK** at residual L2 **6.86e-16**.  Refinement then converges to 166.557530 mK / 218.643409 uA at full residual L2 **6.19e-16**.

The historical 143.777852-uA number is a *transient waveform baseline from a reused known serial steady restart*, not an archived residual-qualified full nonlinear steady solution.  Its excellent COMSOL waveform comparison is real, but it does not establish a 143-uA exact fixed point.

## Answers required by the audit

1. **Circuit equations:** same nominal steady equation and same constants; current's native implementation has different ownership/timing of the update, but no demonstrated algebraic law difference.
2. **TES law:** same affine `R(T,I)` and same clipping.  The selected-point table is `tes_resistance_comparison.csv`.
3. **Joule injection:** same `q=P_relaxed/V_TES`; its raw target is `P_raw=I^2R`. `V_TES=4e-14 m3` in both SIFs and is independently reproduced from both meshes. At the refined current state, `P_raw-P_relaxed=5.17e-15 W`; at Gate3 it is `4.73e-14 W`. The Joule table records this rather than silently treating it as equality.
4. **Thermal path:** not yet proven equivalent.  Materials, nominal dimensions, bath condition, mortar enablement, and all non-Stycast body volumes agree.  Stycast mesh volumes differ by about 0.51% (faceted-cylinder discretization); that alone is not evidence for a 143-to-218-uA shift.  No body/interface flux archive exists for the historical run, so a conductance comparison cannot honestly be claimed.
5. **Stycast COMSOL equivalence:** not established; both use 20 um nominal thickness and identical material properties, but no COMSOL interface-flux or contact-conductance record is present.
6. **Dominant evidence-backed candidate:** the historical 143-uA baseline was accepted/reused without residual qualification, while the current 143-uA checkpoint demonstrably is not a thermal fixed point.  Confidence: **high** for this statement; **low** for assigning the remaining difference to physical thermal conductance.
7. **Code/config to change now:** none in production.  Do not tune HYPRE/GPU and do not alter physics.
8. **Will a correction restore 143 uA?** Unknown.  First replay the historical mesh/state under a residual-capturing MUMPS solve.  Only a residual-qualified 143-uA result justifies a controlled transplant of its circuit update semantics or thermal mesh/interface.

## Minimal next controlled test

Run `mesh_singlepixel_prod_v2` with CPU MUMPS and the current native HeatSolve hook from the exact saved/reused historical restart, capture the full constrained residual plus body/interface fluxes.  This is a diagnostic mesh/coupling transplant, not a production-physics change.  It separates the currently proven non-fixed-point issue from a possible mesh/interface effect.

HYPRE/GPU tuning remains **NO-GO**.
