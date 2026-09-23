# Target-side mesh-density hypothesis: formally rejected

The actual target is `abs` body 100 / boundary 1004, not `SiO2_2`. Refining that
contact from 72 faces (mean edge 50.15 um) to 4,349 faces (mean edge 9.86 um)
changed `G_eff` from `2.020558178452e-8` to `2.030892975143e-8 W/K`, i.e. +0.511%
and moved away from the historical `1.878294461631e-8 W/K`. The patch error
improved, but conductance did not. All three CPU-native MUMPS solves completed
with exit 0 and residuals O(1e-15). This candidate is rejected as the dominant
source of the remaining difference.
