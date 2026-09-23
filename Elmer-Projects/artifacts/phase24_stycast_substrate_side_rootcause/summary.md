# Phase24 Stycast substrate-side controlled test

The TES-side h=10 um density probe is held fixed. Only the final Stycast layer at the Stycast/SiO2_2 contact is additionally refined to h=10 um.

- Historical G_eff: `1.878294461631e-08 W/K`
- TES-side-only baseline G_eff: `2.029939858208e-08 W/K` (ratio `1.080735689`)
- TES+substrate-side G_eff: `2.020558178452e-08 W/K` (ratio `1.075740902`)
- substrate-side improvement: `9.381679756241e-11 W/K`; residual-8% contribution: `0.061866`
- classification: `C`
- interface face area is preserved in the accepted paired mesh; the substrate mesh is unchanged.
- substrate patch errors are in `substrate_patch_tests.csv`; local face stiffness and mortar support are in `substrate_interface_operator.csv`, `mortar_operator_summary.csv`, and `substrate_constraint_support.csv`.
- physics/material/TES law/circuit/bath/geometry dimensions unchanged; production mesh not overwritten.
- full nonlinear run: not performed because the 2–3% historical-equivalence condition was not met.
- HYPRE/GPU: NO-GO.
