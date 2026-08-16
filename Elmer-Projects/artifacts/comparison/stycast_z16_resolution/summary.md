# Stycast vertical-resolution and time-step validation

All crossing levels use the COMSOL peak current drop (7.774844 uA).

| Series | Baseline [uA] | t10 [us] | t50 [us] | Drop at 40 us [uA] | Drop at 100 us [uA] |
|---|---:|---:|---:|---:|---:|
| COMSOL | 143.055049 | 41.7428 | 92.8374 | 0.681418 | 4.247100 |
| Elmer z1, 10 us | 147.375245 | 6.3708 | 43.9263 | 3.649897 | 6.228768 |
| Elmer z8, 10 us | 147.355446 | 36.7784 | 95.3311 | 0.919637 | 4.110400 |
| Elmer z16, 10 us | 147.355947 | 37.8497 | 95.9654 | 0.870323 | 4.082208 |
| Elmer z16, 5 us | 147.355945 | 38.7863 | 93.6845 | 0.839683 | 4.196130 |
| Elmer z16, 2.5 us | 147.355946 | 39.8251 | 92.0008 | 0.786408 | 4.286929 |
| Elmer z16, 1.25 us | 147.355944 | 40.6410 | 92.3019 | 0.743337 | 4.274085 |
| Elmer z32, 1.25 us | 147.355426 | 40.9409 | 92.4487 | 0.727431 | 4.267419 |
| Elmer z32, 0.625 us | 147.355140 | 41.1820 | 92.1898 | 0.714105 | 4.282276 |

## Separation of effects over 0--100 us

- Spatial refinement, z8 to z16 at 10 us: maximum change 0.051696 uA (0.665% of COMSOL peak).
- Time refinement, 10 us to 5 us at z16: maximum change 0.113922 uA (1.465% of COMSOL peak).
- Time refinement, 5 us to 2.5 us at z16: maximum change 0.090800 uA (1.168% of COMSOL peak).
- Time refinement, 2.5 us to 1.25 us at z16: maximum change 0.044140 uA (0.568% of COMSOL peak).
- Spatial refinement, z16 to z32 at 1.25 us: maximum change 0.016049 uA (0.206% of COMSOL peak).
- Time refinement, 1.25 us to 0.625 us at z32: maximum change 0.015007 uA (0.193% of COMSOL peak).
- z16/10 us versus COMSOL: maximum residual 0.215487 uA (2.772% of COMSOL peak).
- z16/5 us versus COMSOL: maximum residual 0.158265 uA (2.036% of COMSOL peak).
- z16/2.5 us versus COMSOL: maximum residual 0.105074 uA (1.351% of COMSOL peak).
- z16/1.25 us versus COMSOL: maximum residual 0.062276 uA (0.801% of COMSOL peak).
- z32/1.25 us versus COMSOL: maximum residual 0.046625 uA (0.600% of COMSOL peak).
- z32/0.625 us versus COMSOL: maximum residual 0.035782 uA (0.460% of COMSOL peak).

## Interpretation

The z8-to-z16 spatial change and successive time-step changes show which
discretization is limiting the early response.  The remaining error is not by itself evidence that the
0.16 um TES needs through-thickness refinement.  With the model properties,
the TES diffusivity is about 3.15 m2/s and its thickness diffusion scale is
only 8.1e-15 s.  In contrast, the 20 um Stycast diffusivity is about
5.09e-7 m2/s, giving a slab diffusion scale L2/(pi2 alpha) of 79.6 us.  A
single z16 Stycast element has h2/alpha = 3.07 us, making the 2.5 us and
1.25 us grids physically relevant convergence checks.
