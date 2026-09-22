# Phase24 thermal-network localization

## Scope and provenance

Three already-built meshes were compared with the same diagnostic semantics: CPU native `HeatSolve`, direct MUMPS, frozen TES source, `P0 = 3.203004762115138e-10 W`, bath `150 mK`, unchanged materials/TES law/circuit constants, and actual native full-system captures. Production meshes and production settings were not overwritten. HYPRE/GPU settings were not used.

## Three-way fixed-power result

| metric | historical | original Phase24 | Phase24 TES–membrane mortar |
|---|---:|---:|---:|
| T at 0.95P (mK) | 166.200114399 | 160.967076184 | 160.967079223 |
| T at 1.00P (mK) | 167.052750880 | 161.544289413 | 161.544292463 |
| T at 1.05P (mK) | 167.905387362 | 162.121502641 | 162.121505702 |
| symmetric dT (mK) | 1.705272963 | 1.154426457 | 1.154426478 |
| G_eff = dP/dT (W/K) | 1.878294462e-08 | 2.774542063e-08 | 2.774542010e-08 |
| G bath-flux derivative (W/K) | 1.878294727e-08 | 2.775528465e-08 | 2.774542614e-08 |
| G_secant at P0 (W/K) | 1.878292121e-08 | 2.774536091e-08 | 2.774535358e-08 |
| bath flux at P0 (W) | 3.203008381e-10 | 3.204150713e-10 | 3.203013019e-10 |

The requested derivative is `dP/dT`; bath flux is reported separately because its tiny native observer mismatch is not the imposed-power definition.

## Localization

- Original and TES–membrane-mortar Phase24 are indistinguishable at this resolution: `ΔG/G = -1.890e-08` and center-temperature difference is `3.049992699e-06 mK`. TES–membrane topology contribution is therefore effectively **0%**.
- Historical total thermal resistance is `5.323978856e+07 K/W`; original Phase24 is `3.604198377e+07 K/W` (`Phase24/historical = 0.676975`). The excess conductance is outside TES–membrane.
- Native path partition is in `native_path_heatflow.csv`. Historical heat is carried by explicit mortar reactions, while original Phase24 has approximately zero TES–membrane mortar reaction and nearly all bath flux is classified as the remaining conformal/shared-node path. The transplanted variant moves that heat into the TES–membrane mortar reaction but leaves T(P) unchanged.
- Bath Dirichlet area is equal: `1.751000000e-05 m2` in all three cases. Assignment is SiO2_2, boundary 30 in historical and 1804 in both Phase24 meshes.
- TES→Stycast and Stycast→substrate physical areas are nearly the same between original and modified Phase24; their native mortar row counts are unchanged. No new direct TES→Stycast, TES→bath-connected, membrane→bath-connected, or Stycast→bath-connected shared-node path was found in the targeted audit.
- `unintended_connection_audit.csv` reports no same-boundary duplicate faces and no interface represented by both shared nodes and native mortar constraints. Internal faces appearing under two body boundary records are treated as expected FE interface bookkeeping.

## Controlled-topology conclusion

The only completed transplant is TES–membrane. It does not move Phase24 toward historical behavior, so the next minimum controlled fix is **membrane→Stycast only**, followed by the same three frozen-power points. Do not change materials, conductivity, circuit constants, production mesh, HYPRE, or GPU settings until that test is complete.

The fixed-power network difference explains the direction of the 143→218 µA branch shift and is a strong root-cause candidate, but it is not by itself a full nonlinear current proof; the full branch must be re-run after the responsible outer interface is isolated.

Full ElmerSolver status: all six new historical/original runs completed with native capture; three modified runs were reused from the prior native MUMPS diagnostic and matched. Physics changes: none. HYPRE/GPU: **NO-GO**.
