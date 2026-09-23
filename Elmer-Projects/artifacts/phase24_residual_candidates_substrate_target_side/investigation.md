# Phase24 residual candidate: actual substrate-side target contact

## Scope

Only the target-side mesh around the actual Stycast-to-substrate contact was changed. The actual contact parent is body `abs` (body 100), boundary 1004; `SiO2_2` is downstream and is not the contact face. TES-side and Stycast-side h=10 um controls, materials, geometry, circuit, power points, CPU-native MUMPS, and HYPRE/GPU NO-GO status were kept unchanged.

The diagnostic target sizing is `SUBSTRATE_CONTACT_TARGET_REFINE_H=10 um`, applied to the `abs` contact footprint and the first target-side 20 um region.

## Mesh evidence

| case | target contact faces | target edge mean / median | target adjacent depth mean | target patch rows |
|---|---:|---:|---:|---:|
| TES+Stycast-side h=10 um control | 72 | 50.15 / 44.22 um | 12.48 um | 2,370 |
| actual target-side h=10 um candidate | 4,349 | 9.86 / 10.00 um | 2.18 um | 2,369 |
| historical reference | 3,333 | 10.50 / 10.37 um | 3.05 um | 1,671 |

The target-side contact is therefore genuinely refined; this is not the earlier ineffective `SiO2_2`-body probe.

## Fixed-power result

| case | G_eff (W/K) | G_secant (W/K) | ratio to historical |
|---|---:|---:|---:|
| historical native CPU MUMPS | 1.878294461631e-8 | — | 1.000000 |
| TES+Stycast-side h=10 um control | 2.020558178452e-8 | — | 1.075741 |
| actual target-side h=10 um candidate | 2.030892975143e-8 | 2.030885859708e-8 | 1.081243 |

The candidate increases G_eff by `1.03348e-10 W/K` versus the control (`+0.511%`) and explains `-0.63%` of the baseline-to-historical gap. It does not reduce the residual; reject as the primary residual candidate.

TES temperatures for the candidate were 0.164982895620 K, 0.165771466165 K, and 0.166560036710 K at 0.95P, 1.00P, and 1.05P. All three CPU-native MUMPS runs finished with `ALL DONE / exit 0`; full residuals were `6.47e-16`, `6.57e-16`, and `6.35e-16`, with constraint residuals below `5e-21`.

## Operator check

For the Stycast-to-substrate mortar rows, the target-side candidate changed the mean support from 10.67 to 18.86 nodes and the mean x/y support spans from 177.8/141.4 um to 31.3/30.9 um. Constant and linear patch errors remained at O(1e-17); the radial-quadratic L2 error improved from `1.997e-9` to `6.968e-12`. Despite this substantial operator improvement, G_eff moved in the wrong direction, so the remaining discrepancy is not explained by this target-side trace discretization alone.

## Decision

Reject this candidate as the dominant source of the residual. Retain the generated diagnostic mesh and fixed-power captures as evidence. The next investigation should move to a mechanism that can leave the target-side patch tests clean while still shifting the network conductance.
