# Phase24 Stycast 32-layer graph diagnostic

- Source mesh: work\meshes\mesh_phase24_stycast_density_10um_plus_substrate_10um
- Diagnostic mesh: work\meshes\mesh_phase24_stycast_density_10um_plus_substrate_10um_layer_graph_diag
- Stycast elements relabelled: 58,453
- Stycast nodes: 13,454
- z range: 0.00019216 .. 0.00021216 m
- Nominal layer thickness: 6.249999999999999e-07 m
- Element-plane crossings: 0 (nonzero means existing elements are not layer-resolved)
- All 32 layers populated: True
- Shared-node edges: 31 total, 0 non-neighbor skip edges
- Element-connectivity edges: 31 total, 0 non-neighbor skip edges
- Cross-body shared-node rows: 0
- Stycast-parent external faces: 9162
- Unexpected lateral/other faces: 0

## Interpretation

The layer graph is a clean series graph under this diagnostic relabelling.

The copied diagnostic mesh is topology-only; no material or solver input was generated.

Source mesh.names was preserved in memory only for name discovery: True
