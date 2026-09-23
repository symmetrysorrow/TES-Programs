# Phase24 Stycast geometry and material integrals

- Stycast conductivity: 2.69094e-06 W/(m K)
- Stycast density: 2400 kg/m3
- Stycast heat capacity: 0.00122 J/(kg K)

| quantity | historical | Phase24 | Phase24 / historical |
|---|---:|---:|---:|
| total volume (m3) | 3.894023393625e-12 | 3.874103524109e-12 | 0.994884502 |
| V / thickness effective area (m2) | 1.947011696813e-07 | 1.937051762054e-07 | 0.994884502 |
| integral k dV (W m) | 1.047858331084e-17 | 1.042498013717e-17 | 0.994884502 |
| integral rho cp dV (J/K) | 1.140170049653e-11 | 1.134337511859e-11 | 0.994884502 |
| 1-D series R (K/W) | 3.817309918089e+07 | 3.836941065242e+07 | 1.005142665 |

## Layer uniformity

- Historical layer effective-area min/max ratio: 1.000000000000
- Phase24 layer effective-area min/max ratio: 0.994684862444
- Historical element-plane crossings: 0
- Phase24 element-plane crossings: 0

The Stycast material is constant-k in the SIF, so the material conductivity integral is exactly k times mesh volume. The volume and effective-area ratios are approximately 0.995, far from the approximately 1.48 conductance excess.
