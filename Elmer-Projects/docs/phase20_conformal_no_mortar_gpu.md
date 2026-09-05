# Phase20 conformal shared-node route (GPU/HYPRE)

## Status

This is an independent validation route for the single-pixel TES--Stycast--absorber stack. It does not replace the existing Mortar + saddle-point/Block-Schur route. The production recommendation remains **PARTIAL** until transient physics and field-level parity are demonstrated.

The implementation is selected by the generated project flag
`elmer_overrides.conformal_shared_node_interfaces = true`. The existing
`conformal_mortar_interfaces` path remains available and its `_MORTAR_PAIRS`
contract is regression-tested unchanged.

## What was implemented

- `generate_project_geometry.py` now has a separate shared-node interface
  retagging path. The three contacts are published as semantic boundary
  pairs:
  `Membrane_SiNx__zmax/TES__zmin`, `TES__zmax/Stycast__zmin`, and
  `Stycast__zmax/abs__zmin`.
- The contact footprints are imprinted/split before extrusion and the final
  Gmsh mesh is converted with `ElmerGrid ... -merge 1e-10`. The route does not
  emit Mortar BCs, Galerkin projection, Block mode, or Schur options.
- `scripts/analysis/check_conformal_interfaces.py` checks the post-ElmerGrid
  mesh itself: shared Elmer node IDs, identical surface-element partitions,
  coordinate gaps, body adjacency, tetrahedral volume validity, and duplicate
  connectivity.
- `scripts/analysis/summarize_conformal_runs.py` records direct, CPU-HYPRE,
  GPU-HYPRE, and Mortar-reference solver evidence. Vectors from different
  meshes are explicitly marked non-comparable rather than being interpolated.

## Mesh/topology gate

The generated mesh is 3D all-tetra and contains:

- 26,705 nodes
- 132,458 tetrahedra
- minimum absolute tetra volume: `1.6092208332073137e-20 m^3`
- zero-volume elements: `0`
- non-finite-volume elements: `0`
- duplicate volume connectivity: `0`

All three interfaces pass the post-conversion gate:

| interface | shared nodes | surface elements | coordinate gap |
| --- | ---: | ---: | ---: |
| Membrane--TES | 79 | 132 / 132, matched 132 | 0 m |
| TES--Stycast | 51 | 81 / 81, matched 81 | 0 m |
| Stycast--absorber | 31 | 50 / 50, matched 50 | 0 m |

Evidence: [interface_connectivity.json](../artifacts/phase20_conformal/interface_connectivity.json).

## Solver gates

### Direct conformal and Mortar reference

The stock Windows Elmer 26.1 executable does not contain MUMPS, so the direct
validation uses UMFPACK. Both runs completed with `ALL DONE`. The Mortar
reference uses a separate `mesh_refined_3x` mesh and a relaxed `1e-6`
reference tolerance; its result vector therefore has a different length and
is not compared entry-by-entry with the conformal mesh.

The direct run is useful as a runtime/syntax gate, but it is not yet a proof
of physical parity: the current Elmer `SaveResult` vector and the UDF
iteration-series scalar values do not expose the same observable in this
setup, and a robust reference iteration-series comparison is still missing.

### HYPRE CPU/GPU

The CPU and GPU cases use FlexGMRES + BoomerAMG with 1,000 maximum iterations
and `1e-7` convergence tolerance. The earlier `1e-10` gate reached a residual
of about `1.94e-8` at the iteration limit and was rejected for this smoke
case. The successful runs used the Phase20 WSL HYPRE-enabled Elmer binary.

Both runs completed with three nonlinear iterations and `ALL DONE`:

- CPU: HYPRE solution-time sum `0.269890054 s`, solver wall `3.57 s`, final
  Result Norm `0.15064936628904274`.
- GPU: HYPRE solution-time sum `0.569229011 s`, solver wall `4.91 s`, final
  Result Norm `0.15064935955824685`; the log records HYPRE IJ matrix/vector
  migration to device memory.
- Final scalar Result Norm difference: about `7.3e-9`.
- Saved full result-vector parity: relative L2 `9.75e-7`, maximum absolute
  entry difference `5.97e-6`.

The GPU is slower on this small smoke mesh, so this is a correctness/integration
gate, not a production performance claim. Evidence:
[run_summary.json](../artifacts/phase20_conformal/run_summary.json).

## Relationship to Block-Schur

The conformal route removes the saddle-point interface unknowns from this
case, so it has no Mortar projector or Schur complement to precondition. That
is structurally simpler than the current Mortar + Block-Schur path. However,
the present evidence does not include an apples-to-apples production-size
benchmark against the Block-Schur artifacts, so no speedup claim is made.

## Remaining gates and recommendation

The route is currently **PASS** for mesh topology and CPU/GPU HYPRE solver
integration, but **PARTIAL / NOT PROVEN** for physical parity. The next gates
are:

1. add a result-field reader for temperature jumps and heat-flux balance at
   all three interfaces;
2. compare the conformal and Mortar references on a common mesh or with a
   documented interpolation/error metric;
3. run a transient one-step and the production-relevant pulse/time grid;
4. repeat the same checks for the dual-TES geometry and larger meshes.

Until those gates pass, keep Mortar + Block-Schur as the production fallback
and promote this route only to the next transient physical-validation stage.

The working tree already contained unrelated Phase20 changes, so this work
was kept on the current `phase20-native-probe-integration` checkout without
resetting or switching branches.
