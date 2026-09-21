# No-mortar / mortar COMSOL parity diagnostic

## Scope

The existing cases are not an identical-mesh A/B. The no-mortar case uses the fine all-tetra conformal mesh; the mortar reference uses the refined-3x all-tetra mesh. This report therefore identifies established differences but does not attribute every difference to mortar alone.

## Mesh and interface topology

| quantity | no-mortar | mortar reference |
|---|---:|---:|
| nodes | 90,872 | 63,172 |
| supported elements | 459,683 | 312,997 |
| element mix | tet:459,683 | tet:312,997 |

Interface face counts (shared finite-element faces; mortar contact faces are explicit boundary pairs and therefore do not appear as internal faces here):

| interface | no-mortar faces | mortar reference faces |
|---|---:|---:|
| Membrane_SiNx / Membrane_Si1 | 690 | 4,642 |
| Membrane_SiNx / SiNx | 136 | 360 |
| SiO2_1 / Membrane_Si1 | 582 | 4,502 |
| SiO2_1 / Si_1 | 20,314 | 13,988 |
| SiO2_1 / Si_2 | 20,312 | 13,974 |
| Si_1 / Membrane_Si1 | 132 | 564 |
| Si_1 / SiNx | 20,430 | 14,032 |
| Si_2 / SiO2_2 | 20,324 | 5,934 |
| TES / Membrane_SiNx | 360 | 11 |
| TES / Stycast | 241 | 0 |
| abs / Stycast | 208 | 0 |

## SIF parity

The common circuit/material/pulse constants are compared in `diagnostic.json`. `Apply Mortar BCs` is intentionally different; the report also records the pulse discrete norm, which is mesh-dependent even though the source is normalized to the requested pulse energy.

## Waveform evidence

| metric | COMSOL | no-mortar | mortar reference |
|---|---:|---:|---:|
| baseline current [µA] | 143.055049 | 154.851455 | 148.130158 |
| minimum current [µA] | 135.280206 | 148.585783 | 140.404956 |
| peak current drop [µA] | 7.774844 | 6.265672 | 7.725202 |
| peak delay [ms] | 0.428000 | 0.100001 | 0.500001 |
| 10–90 rise [ms] | 0.161440 | 0.074037 | 0.174026 |

## Same-mesh mortar probe

| quantity | no-mortar coarse | no-mortar fine | mortar enabled on the same fine mesh | fine→same-mesh mortar |
|---|---:|---:|---:|---:|
| final TES temperature [K] | 0.167254722831 | 0.168249101998 | 0.168249529718 | +4.28e-07 |
| final raw current [µA] | 191.673156340 | 154.858777131 | 154.843477687 | -0.015299445 |
| final relaxed power [W] | 3.91170078802e-10 | 3.3824942464e-10 | 3.38251590184e-10 | +2.17e-15 |

The same-mesh probe converged with MUMPS and changed the raw current by only about -0.0153 µA. The existing coarse→fine no-mortar steady change is much larger (about -36.8 µA), so the dominant observed sensitivity is spatial discretization. The same-mesh mortar flag alone does not reproduce the large improvement seen in the separate mortar reference.

## Current conclusion

The no-mortar waveform remains substantially different from COMSOL in the existing data. The same-mesh probe shows that simply adding mortar constraints does not fix it. The strongest current hypothesis is insufficient spatial resolution/topology in the thin TES/Stycast/interface stack: the no-mortar fine mesh has only 1,615 TES tetrahedra and 901 Stycast tetrahedra, while the mortar refined-3x reference has 7,929 and 2,865 respectively. The next test should refine the no-mortar interface stack itself, then compare the converged MUMPS waveform.

This model solves a thermal PDE plus a lumped TES circuit; it does not solve an electric-potential/current-density PDE. Temperature/heat-flux continuity requires VTU field output and was not assessed here because the stored long runs do not include VTU fields. Electric potential/current-density continuity is not an available field in this model.
