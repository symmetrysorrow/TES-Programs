# Production-v2 no-mortar MUMPS reference

- Case: `case_tes_steady_prod_v2_nomortar`
- Mesh: `mesh_singlepixel_prod_v2`
- Interface: no-mortar
- Initial state: `T_0`; restart: none
- Solver: Stage 11 MUMPS
- Exit: `0`; `ALL DONE`: recorded
- Steady current: `143.53734493231093 µA`
- COMSOL steady current: `143.055049 µA`
- Difference to COMSOL: `0.3371400979429485%`

This is the frozen MUMPS reference for the rewritten no-mortar Gate 3. The
historical mortar reference remains in `prod_v2_mumps_vs_comsol`.
