# Phase24 Stycast geometry and material integrals

- Stycast conductivity: 2.69094e-06 W/(m K)
- Stycast density: 2400 kg/m3
- Stycast heat capacity: 0.00122 J/(kg K)

| quantity | historical | Phase24 | Phase24 / historical |
|---|---:|---:|---:|
| total volume (m3) | 3.894023393625e-12 | 3.876101874674e-12 | 0.995397686 |
| V / thickness effective area (m2) | 1.947011696813e-07 | 1.938050937337e-07 | 0.995397686 |
| integral k dV (W m) | 1.047858331084e-17 | 1.043035757863e-17 | 0.995397686 |
| integral rho cp dV (J/K) | 1.140170049653e-11 | 1.134922628904e-11 | 0.995397686 |
| 1-D series R (K/W) | 3.817309918089e+07 | 3.834964467116e+07 | 1.004624867 |

## Layer uniformity

- Historical layer effective-area min/max ratio: 1.000000000000
- Phase24 layer effective-area min/max ratio: 0.994902990462
- Historical element-plane crossings: 0
- Phase24 element-plane crossings: 0

The Stycast material is constant-k in the SIF, so the material conductivity integral is exactly k times mesh volume. The volume and effective-area ratios are approximately 0.995, far from the approximately 1.48 conductance excess.
