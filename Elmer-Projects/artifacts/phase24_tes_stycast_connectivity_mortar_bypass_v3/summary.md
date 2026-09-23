# Phase24 TES↔Stycast single-interface conforming bypass

## Controlled change

The source is the validated 10 um TES-side Stycast refinement. Geometry resolution, materials, TES law, circuit constants, bath, power points, solver, and nonlinear settings are unchanged. Only TES↔Stycast (`1105/1204`) was changed from mortar to shared-node continuity. Membrane↔TES (`2305/1104`) and Stycast↔abs (`1205/1004`) remain mortar interfaces; only the former contributes active mortar constraints in this solve.

The corrected mesh has 97,076 nodes and 490,831 elements. Boundary strip post-processing removed 727 facets from 1105/1204 only. The post-strip audit retains 490 facets on 1004 and 74 on 1205.

## Result

Using the same symmetric frozen-power slope, `G_eff = (P+ - P-) / (T+ - T-)`:

| case | 0.95P TES [K] | 1.00P TES [K] | 1.05P TES [K] | G_eff [W/K] | historical ratio |
|---|---:|---:|---:|---:|---:|
| historical | 1.662001143989e-1 | 1.670527508805e-1 | 1.679053873620e-1 | 1.878294461631e-8 | 1.000000 |
| refined mortar baseline | 1.649899176177e-1 | 1.657788584203e-1 | 1.665677992230e-1 | 2.029939858208e-8 | 1.080736 |
| TES→Stycast conforming | 1.622559386442e-1 | 1.628548102859e-1 | 1.634500318620e-1 | 2.682374135030e-8 | 1.428090 |

The TES↔Stycast bypass changes `G_eff` by **+32.14%** from the refined-mortar baseline and leaves it **+42.81%** above historical. The center-point `G_secant` is `2.491677971807e-8 W/K`; the large offset shift is secondary to the slope result. This is a strong **NO-GO** for parity and identifies the TES↔Stycast mortar/topology realization as a dominant diagnostic candidate, not as a production fix.

## Actual graph audit

`abs/Pb -> Stycast -> TES -> Membrane_SiNx -> substrate network -> SiO2_2 -> bath`.

- Membrane_SiNx↔TES: mortar, boundary IDs 2305/1104, shared-node count 49, active mortar constraints.
- TES↔Stycast: shared-node, 139 shared nodes and 241 conforming faces, boundary IDs 1105/1204 stripped from the solver mesh.
- Stycast↔abs: mortar, boundary IDs 1205/1004, shared-node count 0.
- Material signatures match the source project; no material or physics parameter was changed.

## Solver and capture

All three corrected runs finished `MAIN: *** Elmer Solver: ALL DONE ***`, exit 0, CPU-native HeatSolve with direct MUMPS. The center capture reports 97,076 primal rows, 203 constraint rows, and 1,334,953 matrix nonzeros. The final center nonlinear relative change is `4.6240616764e-9`.

HYPRE/GPU and full nonlinear production parity remain **NO-GO**. The next useful diagnostic is a controlled physical-network transplant or trace/flux comparison focused on the TES↔Stycast region; conforming TES↔Stycast should not be adopted as the production interface treatment.
