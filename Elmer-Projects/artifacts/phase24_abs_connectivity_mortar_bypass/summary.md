# Phase24 actual `abs` connectivity and single-interface mortar bypass

## Result

The actual thermal graph is **abs/Pb -> Stycast -> TES -> Membrane_SiNx -> the shared membrane/substrate network -> SiO2_2 -> bath**. `abs` is not a bath-facing SiO2 body: it is a 1 mm x 1 mm x 0.7 mm Pb/absorber body above Stycast, connected to Stycast through the mortar pair (`Stycast__zmax` / `abs__zmin`, Phase24 boundary 1004 on the abs side). Historical has the same physical abs body as body 10; Phase24 has it as body 100. Material IDs and parameter values are identical.

The controlled bypass changed **Stycast -> abs only** from mortar to shared-node continuity. The Membrane/TES and TES/Stycast mortar pairs remained active. The bypass used 99,630 nodes and 497,795 elements; the baseline used 99,185 / 495,944. Outer geometry, materials, thickness, power and bath were unchanged.

| case | 0.95P TES [K] | 1.00P TES [K] | 1.05P TES [K] | G_eff [W/K] | G_secant [W/K] | historical ratio |
|---|---:|---:|---:|---:|---:|---:|
| historical | 1.662001143989e-01 | 1.670527508805e-01 | 1.679053873620e-01 | 1.878294461631e-08 | 1.878292121058e-08 | 1.000000 |
| refined mortar baseline | 1.649899176177e-01 | 1.657788584203e-01 | 1.665677992230e-01 | 2.029939858208e-08 | 2.029934407666e-08 | 1.080736 |
| Stycast->abs conforming | 1.665505837435e-01 | 1.673409941716e-01 | 1.681245645209e-01 | 2.034970698654e-08 | 1.847071010133e-08 | 1.083414 |

`f_mortar = (G_base-G_conf)/(G_base-G_hist) = -0.033175`. On the requested differential `G_eff` metric, the bypass changes the result by +0.25% and leaves the case at +8.34% relative to historical. Thus this is **Case B**: Stycast/abs conforming continuity does not remove the remaining ~8% differential conductance gap. The center-point `G_secant` shifts substantially, so the bypass also exposes an absolute thermal-offset/network effect, but it is not evidence that the mortar operator alone explains the slope mismatch.

## Required decisions

1. `abs` role: Pb absorber, body 10 historical / body 100 Phase24, volume 7.0e-10 m3, bbox [-0.5,0.5] mm x [0.5,1.5] mm x [0.21216,0.91216] mm, thickness 0.7 mm, 1 mm2 footprint.
2. Historical corresponding body/path: yes; same Pb material and geometry. Historical uses fully mortar-separated TES/Stycast/abs contact boundaries; Phase24 has the same named bodies but a different mesh/contact surface realization and has shared nodes on the Membrane/TES contact.
3. Material sequence: same named materials and SIF values. The actual route is not a direct `TES -> Stycast -> abs -> SiO2_2` vertical stack; the bath path continues from TES through the membrane/substrate network to SiO2_2.
4. Bypass target: Stycast <-> abs only. Mesh resolution was not the manipulated variable.
5. Trace-to-flux: historical and refined-mortar low-mode patch tests remain clean for constant/linear modes and show the known radial-quadratic error; the bypass has exact shared-node trace continuity by construction, but a separate imposed-mode reaction capture was not run. See `trace_to_flux_modes.csv`.
6. Segment resistance: existing trace/body decomposition and the bypass body-average fallback are in `segment_resistance.csv`. The bypass shared Stycast/abs trace has no mortar jump; its local flux distribution was not separately captured.
7. Target-side mesh refinement: formally rejected; see `target_side_refinement_rejection.md`.
8. Remaining dominant candidate: coupled physical network realization/topology and the interaction of the remaining TES/membrane and TES/Stycast mortar operators with historical boundary roles. The single abs bypass proves mortar is a contributor but not the sole explanation.
9. Full nonlinear steady: **NO-GO** until a controlled transplant or a second single-interface test brings the frozen-power conductance within 2--3% without changing physics.
10. HYPRE/GPU: **NO-GO**.

Solver status: all three corrected bypass runs finished `ALL DONE`, exit 0, CPU-native HeatSolve with direct MUMPS. The capture reported 99,630 primal rows and 2,524 remaining constraint rows. Physics/material/TES/circuit/bath/production mesh changes: none.

Curated files are listed in the task workspace; large raw matrix captures and generated mesh files are intentionally not part of the curated commit.
