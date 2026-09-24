# Phase24 native branch heat-flow and trace-to-flux operator audit

This audit is diagnostic only. No mesh, SIF, material, TES law, circuit constant, bath BC or production file was changed, and ElmerSolver was not re-run. The audit re-reads the existing opt-in native captures (`Phase24 Full Restriction Capture`). Each capture is the exact CPU/MUMPS saddle-point system `[[K, C^T],[C, 0]][T; λ] = [f; 0]` and the solution MUMPS returned. The captures cover 0.95/1.00/1.05 P0 for three meshes:

| case | mesh | captures | G_eff (W/K) | / historical |
|---|---|---|---:|---:|
| historical | `mesh_singlepixel_prod_v2` | `phase24_thermal_network_localization/capture/historical_*` | 1.878294462e-08 | 1.000000 |
| best Phase24 (Membrane/substrate h=10 µm control) | `mesh_phase24_trace_membrane_substrate_historical` | `p24trace/ms/capture/ms_hist_*` | 1.959750658e-08 | **1.043367** |
| refined mortar parent | `mesh_phase24_stycast_density_10um` | `p24d10b/capture/phase24_historical_density_*` | 2.029939858e-08 | 1.080736 |

G_eff is recomputed from the captures with the harness definition (symmetric 0.95/1.05 P, TES element-node mean). It reproduces the reported historical value exactly. The best value is 1.9e-5 above the earlier 1.959713948e-8, and the ratio 1.043367 is close to the earlier 1.043348. The frozen-power system is **linear**: the matrix A is bit-identical at the three power points, with k(T) evaluated at the 0.16857 K initial condition. Every symmetric differential quantity therefore equals its secant value.

## Method (native, no proxies)

* **Native checks.** The primal residual is 1e-16. The bath flow is the native edge flow into the Dirichlet rows. It closes against the source to 1.1e-6 (historical), 3.3e-5 (best) and 3.3e-6 (parent). These closure errors are the MUMPS balance precision: T carries a 0.15 K offset while ΔT is about 0.017 K.
* **TES→Membrane flow.** Native edge flow across a closed cut around TES+Stycast+abs, plus the mortar reaction `C^T λ`, split into mortar and conformal parts.
* **Branch flow.** For each suspended sheet, the element reaction `r = K_sheet (T − T_bath)`. The element matrices are recomputed and must reproduce the native K entry by entry: 3e-13 relative for historical, and ≤1.3e-6 for Phase24, where the differences sit only in cold frame rows. The reaction is partitioned into TES inflow, inter-sheet exchange and perimeter exit. The per-sheet balance residual is ≤1e-6 P.
* **Branch conductance.** The dissipation conductance `G_i = Σ_e T_e^T K_e T_e / ΔT_src²` is an exact identity, `Σ_i G_i = G_src`. It equals `∂G/∂ln k_i`, the first-order sensitivity of G to branch i's discrete operator. Every element contribution is ≥ 0, so the split has no circulation artefacts.
* **Edge-level layer attribution was rejected.** Splitting lateral cuts by layer produces ±8 P of edge circulation in the thin-sheet tetrahedra (obtuse-tet positive off-diagonals). It is kept in `native_branch_flux.csv` only as `unavailable`.
* **Operators.** They are condensed on sparse local subdomains with sparse LU. No global dense matrix is formed.

## Actual branch graph (`actual_branch_graph.json`)

The previous `Si1/SiNx` versus `Si2/SiO2_2` parallel-branch picture does not match the real adjacency. **The window below the membrane is back-etched.** SiO2_1 window nodes are not shared with Si_2; Si_2's maximum temperature is 0.1500008 K while the SiO2_1 window runs at about 0.1657 K.

```
TES (source; isothermal within 6 µK)
 ├─ Stycast → abs            dead end, net Q ≈ 1e-18 W
 └─ TES/Membrane trace (z=192 µm, TES footprint)
      historical: mortar, 2218 rows; Phase24: conformal (best has 18 extra hanging-node rows)
      └─ three suspended window sheets in parallel, strongly vertically coupled:
           A  Membrane_SiNx  (z191–192, 1 µm)  ──┐
           B  Membrane_Si1   (z176–191, 15 µm) ──┼── window perimeter s = 350 µm → frame (SiNx/Si_1/SiO2_1/Si_2), ≈ bath
           C  SiO2_1 window  (z175–176, 1 µm)  ──┘      frame→bath: 7.5e-4 W/K, i.e. ~1.3e3 K/W, negligible
```

The frame sits at about 0.1500001 K. All thermal resistance lives in the TES footprint plus the suspended annulus between s = 250 and s = 350 µm.

## Native branch heat flow

Q/P is identical at 0.95, 1.00 and 1.05 P0 because the system is linear. The absolute Q values are in `native_branch_flux.csv`.

| branch | quantity | historical | best Phase24 | parent |
|---|---|---:|---:|---:|
| B0 TES→Membrane | Q/P (mortar part) | 1.000001 (1.000) | 0.999967 (0.0343) | 1.000003 (0) |
| B0s Stycast/abs | Q (W) | 8e-18 | −1e-18 | 5e-19 |
| A Membrane_SiNx | TES in / → B / exit | 1.000 / 0.970 / **0.0299** | 1.000 / 0.993 / **0.0073** | 1.000 / 0.970 / 0.0301 |
| B Membrane_Si1 | → C / exit | 0.522 / **0.448** | 0.535 / **0.458** | 0.523 / 0.447 |
| C SiO2_1 window | exit | **0.522** | **0.535** | 0.523 |
| M frame→bath | Q/P | 1.000001 | 0.999967 | 1.000003 |

At P0, Q_exit is A 9.565e-12 / 2.346e-12 W, B 1.4349e-10 / 1.4663e-10 W and C 1.6724e-10 / 1.7131e-10 W (historical / best). Historical mortar injects 35% of P at membrane master DOFs outside the TES footprint, at s = 250–264 µm (`raw/`).

## Branch conductance (`branch_differential_conductance.csv`)

| branch | G_dissipation hist | best | ratio | excess (W/K) | own dQ/dΔT hist → best |
|---|---:|---:|---:|---:|---:|
| A Membrane_SiNx | 9.633e-10 | 1.189e-09 | 1.235 | +2.26e-10 | 5.61e-10 → 1.44e-10 |
| B Membrane_Si1 | 9.298e-09 | 9.461e-09 | 1.018 | +1.63e-10 | 8.51e-09 → 9.05e-09 |
| C SiO2_1 window | 8.493e-09 | 8.943e-09 | 1.053 | **+4.50e-10** | 1.00e-08 → 1.07e-08 |
| TES | 2.65e-11 | 1.55e-12 | 0.058 | −2.50e-11 | – |
| frame + bath network | 2.3e-12 | 2.0e-12 | – | ≈ +2e-13 | – |

The "own" G is the exit flow over the sheet's own ΔT. It is not additive, because the sheets exchange 50–99% of P vertically. Use the dissipation G for accounting.

## Remaining-gap accounting (`remaining_gap_accounting.json`)

The gap is ΔG_eff = +8.146e-10 W/K (+4.337%).

| category | share of gap |
|---|---:|
| TES→Membrane operator (TES body + coupling) | −3.1% |
| branch A Membrane_SiNx | +27.8% |
| branch B Membrane_Si1 | +20.0% |
| **branch C SiO2_1 window** | **+55.2%** |
| merge/frame→bath network | +0.03% |
| **explained** | **99.95%** |
| unexplained | 0.05% (4.1e-13 W/K = G_eff/G_src definition 2.5e-13 + MUMPS balance 6.6e-13) |

In resistance units at fixed P (`zone_resistance_accounting.csv`), best dissipates less than historical in the **TES-edge band s = 230–270 µm**:

* Membrane_Si1 −2.50e6 K/W
* SiO2_1 −0.84e6 K/W
* Membrane_SiNx −0.39e6 K/W

That band totals −3.8e6 K/W, about 170% of the net ΔR = −2.21e6 K/W. It is partly offset by extra dissipation in the outer annulus and at the window perimeter.

## Operator evidence

**Uniform-flux stack NtD, the decisive like-for-like test (`operator_gap_decomposition.json`).** The same network is loaded by a uniform heat flux over the TES footprint (face-sampled, so node snapping does not matter). There is no isothermal-TES edge singularity in this test.

* G_uniform: historical 1.10001e-8, best 1.10907e-8. The ratio is **1.0082**, far below the coupled 1.0434.
* The exact multiplicative split is 1.04337 = 1.00824 (bulk stack operator) × 1.03484 (isothermal TES-edge loading). In log shares that is **19% bulk and 81% TES-edge loading**.
* The parent→best controlled improvement is 92% edge loading.
* Uniform-flux mode energies (historical/best): constant 1.0082, linear-x 1.0198, linear-y 1.0194, radial-quadratic 1.0141.

**Sheet operators in isolation (`membrane_substrate_operator.json`).** Each sheet alone is loaded by uniform footprint heat with an isothermal perimeter, and compared with a 2-D FD reference S = 11.770, converged to 2e-5.

| case | Membrane_SiNx | Membrane_Si1 | SiO2_1 |
|---|---:|---:|---:|
| historical | +0.28% | +0.28% | +0.28% |
| best | **+15.2%** | +0.31% | +0.30% |
| parent | +8.1% | +1.0% | +1.1% |

The B and C in-plane operators in best are as accurate as historical. The stiffness error of Membrane_SiNx (+15%, ×G_A/G ≈ 6%, ≈ 0.9%) accounts for the bulk 0.82%.

**TES/Membrane local condensed operator.** The subdomain is the Membrane_SiNx layer, with the membrane-side trace prescribed and the layer grounded where it meets B, the SiNx frame and Si_1.

| case | trace DOFs | G (W/K) | row sum per area (median) |
|---|---:|---:|---:|
| historical | 2313 | 7.586e-6 | 28.26 |
| best | 180 | 7.687e-6 | 28.26 |

The median row sum per area equals k/t exactly in both meshes. Normalized mode responses (best/historical):

| constant | linear-x | linear-y | radial-quadratic |
|---:|---:|---:|---:|
| 1.014 | 1.012 | 1.020 | 1.019 |

This layer carries only about 0.25% of R, so its operator difference is at most about 0.004% of G. **The historical slave-side (TES) Dirichlet trace map is ill-posed.** C restricted to master DOFs has σ_min/σ_max = 3.7e-11, with 15 near-null directions. Only master-side and NtD operators are used for historical.

**TES-side NtD of the full native saddle A (TES source modes).** Energy ratios (historical/best):

| constant | linear-x | linear-y | radial-quadratic |
|---:|---:|---:|---:|
| 1.0433 (= G ratio) | 0.996 | 0.991 | 1.023 |

The linear modes are dominated by in-plane TES conduction (k = 68) and do not show a coupling defect.

**Network DtN from the membrane-side trace.** It reproduces P from the actual trace: 1.000001 (historical) and 0.99997 (best). The low-frequency generalized spectrum agrees within 0.4%:

| | μ1 | μ2 | μ3 |
|---|---:|---:|---:|
| historical | 0.04031 | 0.10619 | 0.10630 |
| best | 0.04039 | 0.10581 | 0.10588 |

The difference is therefore high-frequency and edge-localized. The constant-mode totals, 1.9757e-8 versus 1.9599e-8, are not like-for-like: the historical trace includes the mortar master ring out to s = 264 µm.

**Mesh at the TES edge band (`tes_edge_band_mesh_statistics.json`).** In best, planes z = 191, 176 and 175 have 10.7–12.7 µm in-plane edges, the same as historical (11.8 µm). The **z = 192 TES/Membrane trace is 38.6–39.0 µm, with 48 nodes on the TES edge line**, against 11.7 µm and 174 nodes in historical. Membrane_SiNx therefore joins a 39 µm top to a 10.7 µm bottom across 1 µm, which also produces the +15% sheet stiffness.

## Conclusions

1. **Downstream branch native Q was obtained directly.** Each sheet's reaction is native-verified with closure ≤ 1e-6 P. TES→Membrane, the dead end and the bath flow use native edge flow plus `C^T λ`. Edge-level layer attribution is unavailable because of circulation.
2. **Largest branch G difference.** In absolute terms it is branch C, the SiO2_1 window sheet (+4.50e-10, ×1.053). In relative terms it is branch A, Membrane_SiNx (×1.235).
3. **The 46.29% controlled improvement is explained at branch level.** Parent→best: B Membrane_Si1 +93%, C +23%, A −17%. Operator split: 92% TES-edge loading.
4. **The TES/Membrane local condensed operator differs by only 1.2–2.0% per mode.** That component carries about 0.25% of R, and the TES body plus coupling share is −3.1%.
5. **Low-order modes.** The largest mesh sensitivity is in the linear and radial flux modes of the stack: 1.4–2.0% against 0.8% for the constant mode. The coupled gap, however, is carried by the edge-concentrated loading.
6. **The old 45.8% TES→Membrane resistance share is not supported as a coupling or operator effect.** It came from averaging the historical membrane trace (bid 23) over the entire membrane top, including the free annulus. The true TES-trace jump is 1.6 µK in historical and 0.06 µK in best. The root cause does sit at the TES edge, but in how the **z = 192 trace resolution** feeds the stack below it.
7. **Direct evidence explains 99.95%** of the +4.337% gap: exact additive dissipation accounting, independently split 19%/81% by the operator NtD.
8. **Unexplained remainder** is 0.05%, or 4.1e-13 W/K.
9. **Next minimal change:** refine the z = 192 TES/Membrane conformal trace in the TES edge band. See `next_controlled_fix.md`.
10. **Parity 2–3% looks reachable.** The prediction is 1.00–1.012 if the edge-loading term is removed, but it has not been run.
11. **Full nonlinear steady: NO-GO.** No diagnostic fix exists yet that satisfies |G/G_hist − 1| ≤ 0.03.
12. **HYPRE/GPU: NO-GO.**

## Caveats

* Branch shares depend on units. In conductance units, C 55% / A 28% / B 20%. In resistance units at fixed P, B 78% / C 36% / A −17%. Both splits are exact, but they answer different questions, because strongly coupled sheets are not independent parallel branches.
* The dissipation split is a first-order sensitivity attribution, not a one-factor remesh. Its causal localization comes from the uniform-flux versus isothermal-edge operator split and the mesh-band statistics.
* Stycast and abs have G_dissipation ≈ 3e-17, about 1.5e-9 of G. For *steady* G parity, changes to the Stycast mesh are therefore irrelevant. The earlier gate that rejected meshes over Stycast trace changes is stricter than steady parity requires, although it still matters for transients.

Reproduce with `py -3.14 scripts/support/analyze_phase24_native_branch_operator_audit.py` followed by `py -3.14 scripts/support/materialize_phase24_native_branch_operator_audit.py`. Tests are in `tests/test_phase24_native_branch_operator_audit.py`.
