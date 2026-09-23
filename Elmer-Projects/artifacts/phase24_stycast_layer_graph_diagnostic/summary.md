# Phase24 Stycast 32-layer graph diagnostic

- Source mesh: work\meshes\mesh_singlepixel_gpu_fine_stycast32_mortar
- Diagnostic mesh: work\meshes\mesh_phase24_stycast32_explicit_layers_diag
- Stycast elements relabelled: 31,885
- Stycast nodes: 6,193
- z range: 0.00019216 .. 0.00021216 m
- Nominal layer thickness: 6.249999999999999e-07 m
- Element-plane crossings: 0 (nonzero means existing elements are not layer-resolved)
- All 32 layers populated: True
- Shared-node edges: 31 total, 0 non-neighbor skip edges
- Element-connectivity edges: 31 total, 0 non-neighbor skip edges
- Cross-body shared-node rows: 0
- Stycast-parent external faces: 315
- Unexpected lateral/other faces: 0

## Interpretation

The layer graph is a clean series graph under this diagnostic relabelling.

The copied diagnostic mesh is topology-only; no material or solver input was generated.

Source mesh.names was preserved in memory only for name discovery: True
