# Phase24 outer thermal path localization

Actual body connectivity is TES→Membrane_SiNx and TES→Stycast at the TES, then the outer Stycast edge is Stycast→SiO2_2. The bath Dirichlet condition is on SiO2_2. The label TES zmax–Stycast zmin means the TES↔Stycast edge; it is not a membrane edge. See physical_interface_map.csv.

| route | G_eff=dP/dT (W/K) | R_total (K/W) | R/historical |
|---|---:|---:|---:|
| historical | 1.878294462e-08 | 5.323978856e+07 | 1.000000 |
| Phase24 | 2.774542063e-08 | 3.604198377e+07 | 0.676975 |
| Stycast→substrate orientation variant | 2.774542204e-08 | 3.604198194e+07 | 0.676975 |

Reversing only the Stycast→substrate mortar orientation changes Phase24 G_eff by 5.095e-08; this is a null sensitivity and does not move Phase24 toward historical. Therefore Case B applies.

segment_heatflow.csv contains native C-lambda reactions where the validated captures provide them, plus a separately labelled fixed-power global-closure row for the conformal/shared-node remainder. Nonconforming one-sided gradients were not used as conserved heat flows. Consequently, a unique per-interface R_th is not claimed where the path is shared-node/conformal; segment_thermal_resistance.csv reports only directly supportable values. body_temperature_comparison.csv contains TES, membrane, Stycast, and substrate temperatures.

The Phase24 Stycast representation has 32 nominal layers but one converted Elmer body ID. The body-level graph shows no direct Stycast→bath shared-node shortcut or extra parallel Stycast body. A literal layer-by-layer proof is not possible from this converted mesh alone; stycast_layer_graph_audit.csv records this limitation. Confidence for absence of an internal layer shortcut is low-to-moderate.

Full ElmerSolver status: all three controlled variant points completed with CPU native HeatSolve, direct MUMPS, frozen power, and full restriction capture. Solver logs report MUMPS and 189 active restriction rows at each point; no solver failure was reported. Physics/materials/TES law/circuit and production mesh: unchanged. HYPRE/GPU: NO-GO.

For the controlled topology map, 139 TES↔Stycast rows are unchanged and the 189 total restriction rows imply 50 Stycast↔substrate rows; this count is recorded as an inference from the native capture row total. The controlled individual Cλ reaction partition was not materialized after the solver run; its heat-flow cell is therefore NaN rather than guessed.

Dominant resistance conclusion: the Stycast→substrate mortar orientation is not the source of the approximately 1.48× conductance. The remaining difference is localized to outer-network geometry/connectivity or the fused 32-layer representation, but the present evidence does not uniquely select internal Stycast versus substrate-side geometry. Minimum next controlled fix: diagnostic-only per-layer Stycast body labeling and graph audit, followed by the same 0.95/1.00/1.05P rerun; do not change production.
