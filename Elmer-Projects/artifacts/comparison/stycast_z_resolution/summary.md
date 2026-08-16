# Stycast through-thickness resolution validation

The established Phase23 hybrid case uses one prism element through the
20 um Stycast thickness.  The validation case uses eight equal layers
(2.5 um each), with the same Phase19 in-plane fields, material properties,
pulse, timestep schedule, BDF1 lumped mass, circuit model, strict nonlinear
tolerance, solver build, and pinned UDF.

The regenerated one-layer control reproduces the established mesh counts
exactly (223,718 volume elements; 1,674 Stycast prisms).  The eight-layer
mesh has 235,263 volume elements and 13,392 Stycast prisms.  All non-absorber
body counts are identical; the absorber differs by 173 tetrahedra (0.137%).

| Crossing level (relative to COMSOL peak) | COMSOL [us] | Elmer z1 [us] | Elmer z8 [us] |
|---|---:|---:|---:|
| 10% | 41.743 | 6.371 | 36.778 |
| 50% | 92.837 | 43.926 | 95.331 |
| 90% | 203.183 | 136.876 | 203.847 |

Over pulse +0--220 us, the maximum baseline-corrected current difference
from COMSOL falls from 2.9896 uA (38.45% of the COMSOL peak) for z1 to
0.2382 uA (3.06%) for z8.

This result confirms that the sharp early Elmer rise is predominantly a
through-thickness discretization artifact in the low-diffusivity Stycast
layer.  Eight layers restore the distributed diffusion delay while leaving
the steady operating point essentially unchanged.
