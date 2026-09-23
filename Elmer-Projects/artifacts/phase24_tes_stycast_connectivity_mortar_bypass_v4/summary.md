# Phase24 refined TES↔Stycast conforming vs mortar comparison

## Controlled comparison

The source is the validated 10 um TES-side Stycast refinement. Materials, geometry recipe, TES law, circuit constants, bath, power points, solver, and nonlinear settings are unchanged. Only TES↔Stycast coupling is changed from mortar to shared-node continuity. Membrane↔TES and Stycast↔abs remain mortar.

Unlike the earlier coarse conforming test, the TES side now has an explicit `TES_INTERFACE_REFINE_H=10 um` refinement. The conforming interface has 4,581 faces and physical area `1.9472989574e-7 m2`; the refined mortar baseline has 4,579 Stycast-side faces and the same physical area. The v4 mesh has 104,361 nodes and 517,655 elements. Boundary strip removes only TES↔Stycast IDs 1105/1204; Stycast↔abs IDs 1205/1004 remain.

## Conductance

`G_eff = (P+ - P-) / (T+ - T-)`:

| case | 0.95P TES [K] | 1.00P TES [K] | 1.05P TES [K] | G_eff [W/K] | historical ratio |
|---|---:|---:|---:|---:|---:|
| historical | 1.662001143989e-1 | 1.670527508805e-1 | 1.679053873620e-1 | 1.878294461631e-8 | 1.000000 |
| refined mortar | 1.649899176177e-1 | 1.657788584203e-1 | 1.665677992230e-1 | 2.029939858208e-8 | 1.080736 |
| refined conforming TES→Stycast | 1.667347818286e-1 | 1.675347134777e-1 | 1.683275182535e-1 | 2.011007416014e-8 | 1.070656 |

The refined conforming result is **−0.93%** from the refined mortar baseline and **+7.07%** from historical. The simple normalized mortar contribution estimate is `0.12485` of the baseline-to-historical gap. The center-point `G_secant` changes by about −10.0%, so the absolute offset/network state still differs even though the differential slope is nearly formulation-independent at matched interface resolution.

## Actual graph and solver

- Membrane_SiNx↔TES: mortar; boundary IDs 2305/1104.
- TES↔Stycast: shared-node; 2,370 shared nodes and 4,581 conforming faces before boundary strip.
- Stycast↔abs: mortar; boundary IDs 1205/1004; shared-node count 0.
- Solver: CPU-native HeatSolve with direct MUMPS; all three runs exit 0 and report `ALL DONE`.
- Center full restriction capture: 104,361 primal rows, 874 remaining constraint rows, 1,462,231 matrix nonzeros. Final nonlinear relative change: `1.1467390447e-9`.

## Decision

The coarse conforming result remains **NO-GO** because it conflated coupling formulation with interface coarsening. The refined conforming method is **GO as a matched-resolution diagnostic comparator**: mortar vs conforming differs by less than 1% in `G_eff`. It is not a production-interface change. The remaining approximately 7% historical gap is not explained by TES↔Stycast coupling alone and should be pursued in the physical network realization/topology and historical-vs-Phase24 body/interface mapping.
