# Phase24 steady-branch diagnosis

## Result

- Gate3 one-shot MUMPS thermal solve: 168.569115 -> 162.369393 mK; field-re-evaluated raw circuit current 143.568115 -> 390.165660 uA.
- The Gate3 assembled power was 3.199023057e-10 W, only 4.026e-15 W below the saved 3.199063315e-10 W. The 6.20-mK direct move is therefore strong linear-defect evidence (Hypothesis A).
- The intended checkpoint circuit state was not reloaded: both diagnostic assembly rows show the T0/143.537-uA state. Consequently the refined-218-uA one-shot is not a valid nonlinear fixed-point test; branch multiplicity and the existence of a 143-uA high-accuracy fixed point remain open.
- No physical, mesh, mortar, HYPRE-preconditioner, GPU, or production-timestep parameter was changed.

## Endpoint and comparison

- The solver schedule reaches 39.376 us post-pulse, while both series stop at 38.751 us: the final series row is missing.
- Recomputed three-way common-grid triangle check: True. Earlier contradictory values used non-identical comparison/baseline pipelines.
