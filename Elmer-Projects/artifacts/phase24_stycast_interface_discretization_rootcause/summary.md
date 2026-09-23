# Phase24 Stycast interface discretization root-cause audit

This diagnostic directory contains read-only measurements and existing native CPU/MUMPS fixed-power evidence. Physics parameters, materials, TES law, circuit constants, and production meshes were not changed.

## Interface density and fixed-power result

| case | Stycast TES-interface faces | mean face edge (m) | median adjacent depth (m) | G_eff (W/K) | ratio |
|---|---:|---:|---:|---:|---:|
| historical | 3214 | 1.064244498e-05 | 3.125000000e-07 | 1.8782944616314638e-08 | 1.0 |
| phase24_original | 241 | 4.095627221e-05 | 1.562500000e-07 | 2.7745420625771133e-08 | 1.4771603277620164 |
| phase24_h5um | 8308 | 4.923111237e-06 | 1.562500000e-07 | 2.0328561515006625e-08 | 1.0822883168887951 |
| mesh_phase24_stycast_density_10um | 4579 | 9.775442017e-06 | 1.562500000e-07 | 2.0299398582083175e-08 | 1.0807356885060164 |
| mesh_phase24_stycast_density_15um | 2053 | 1.450021074e-05 | 1.562500000e-07 |  |  |
| mesh_phase24_stycast_density_20um | 1165 | 1.915602435e-05 | 1.562500000e-07 |  |  |

## Interpretation

- Existing Stycast-only h=5um control preserves the original Phase24 TES element count and TES-side interface tessellation, yet reduces G_eff from 2.77454e-8 to 2.03286e-8 W/K (1.0823 historical ratio).
- This confirms that the dominant sensitivity is on the Stycast contact-side neighborhood, but h=5um changes local bulk element stiffness and mortar trace/projection resolution together.
- The current evidence therefore does not justify choosing bulk stiffness alone or mortar projection alone as the root cause.
- The 10um probe was selected as the closest geometry-density candidate and its three native MUMPS frozen-power cases were run; it plateaus at 1.080736 historical ratio.
- HYPRE/GPU remains NO-GO for this phase.

## Required decisions

1. Provisional historical-equivalent Stycast contact density: target h=10um, actual mean face edge 9.7754um, 4,579 Stycast-side faces, and 2,414 total native constraints (2,369 classified TES-Stycast). It is the closest tested edge-density probe, not an exact topology match.
2. G_eff does not converge to historical: 2.029939858e-8 W/K, ratio 1.0807357. This is Case B; about 8% remains.
3. The solved sequence original -> h=5um -> h=10um shows a clear first reduction then a plateau, but a complete 3-5 point G_eff(h) curve is not established because 15/20um were geometry-only probes.
4. TES-side refinement is unnecessary for the dominant reduction: the Stycast-only h=5um control leaves TES elements/faces unchanged and matches the TES+Stycast refinement within about 0.20%.
5. Bulk near-interface stiffness versus mortar trace resolution is not fully separable in this generator. Both change together; the controlled evidence supports Stycast-side interface-neighborhood discretization as the cause, without selecting one sub-mechanism.
6. Local static-condensed trace stiffness is reported in local_effective_stiffness.csv. Aggregate condensed trace K is 2.0374e-10 W/K original, 8.0176e-9 h=5um, 3.7951e-9 h=10um, versus 1.9389e-8 historical; this is a local diagnostic, not global G_eff.
7. B*u patch tests show floating-point constant/linear residuals; the radial-quadratic residual is 6.75e-10 original versus 3.99e-12 historical, while h=5um/10um are about 6.0e-10. A low-frequency projection difference exists, but it is not by itself a conductance proof.
8. Thermal discretization explains the observed 47.7% original conductance excess only partially: the Stycast-side change removes about 26.8% of G and leaves about 8.1% above historical.
9. Full nonlinear refined MUMPS was not run because the conditional historical-equivalence criterion was not met; nonlinear_refined_steady.json records this explicitly.
10. The 143 -> 218 uA causal chain is not confirmed by this campaign; current remains an open nonlinear/thermal-coupling question after the Case B plateau.
11. No production mesh replacement is authorized. The minimum next diagnostic is a Stycast contact-side local refinement around target h~10um, followed by a controlled nonlinear run; production choice remains deferred.
12. HYPRE/GPU: NO-GO.

## Solver and change status

- Full ElmerSolver: h=10um 0.95/1.00/1.05P all completed with exit 0; center native-after full residual L2=6.5518e-16, primal L2=6.5196e-16, constraint L2=5.1883e-20.
- Physics/material/circuit/TES law/production mesh: unchanged.
- Native center mortar reaction diagnostic: ||B^T lambda||_2=2.6174e-13, TES-weighted=-2.6844e-18 W; D block is exactly zero in the captured saddle system.

## Requested artifact status

`patch_test_results.csv` applies constant, linear-x, linear-y, and radial-quadratic fields directly to the signed captured B block. Constant/linear residuals stay at floating-point scale; the radial low-order mode is materially larger for original Phase24 than historical, while h=5um and 10um are similar.
