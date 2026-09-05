# Phase20 conformal / Mortar physical parity

## Scope

This is the clean branch `phase20-conformal-physical-parity-latest`, based on
current `origin/main` `d04102b4`. The native-probe branch remains untouched.
The existing Mortar + Block-Schur route is the control; this branch validates
the conformal/shared-node primal route.

## Geometry and topology

The smoke conformal mesh has 26,693 nodes and 132,357 tetrahedra. Its three
interfaces pass with shared nodes, matched surface partitions, connected
element adjacency, and zero duplicate/zero-volume elements.

The previously available 90,872-node production candidate was not valid for
the current branch. Its failure was classified as:

- Membrane/TES: missing semantic left surface (`RETAGGING_ERROR`) plus node
  and partition mismatch
- TES/Stycast: node merge, contact partition, surface mesh, and coordinate-gap
  mismatch
- Stycast/absorber: node merge and surface partition mismatch

Regenerating the same production-size recipe with the current generator fixes
the failure:

- 90,872 nodes / 459,683 tetrahedra
- Membrane/TES: 205 shared nodes, 360/360 matched facets
- TES/Stycast: 139 shared nodes, 241/241 matched facets
- Stycast/absorber: 114 shared nodes, 208/208 matched facets
- coordinate gaps: 0
- zero volume, nonfinite volume, duplicate connectivity: 0

The production topology blocker is therefore **CLOSED**. No GPU benchmark has
been started yet.

## Stycast ideal / OCC / FEM convergence

The OCC BREP kernel volume for the Stycast cylinder is
`3.8956377223044e-12 m3`, equal to the ideal 498 um cylinder. The remaining
error is polygonal Gmsh discretization, not an OCC boolean-volume mismatch.
ElmerGrid preserves the pre-Elmer Gmsh volume to below `8e-13` relative error.

| level | nodes | tetrahedra | Gmsh/FEM Stycast volume error vs ideal |
| --- | ---: | ---: | ---: |
| coarse | 90,872 | 459,683 | -1.064920% |
| medium | 225,821 | 1,170,252 | -0.468533% |
| fine | 382,507 | 1,978,314 | -0.266416% |

The three-level geometry convergence gate is **PASS**: the error decreases
monotonically and Elmer FEM preserves the Gmsh volume. This closes the earlier
“unknown Stycast geometry difference” blocker, although the ideal-cylinder
error is still reported rather than hidden.

## Electrical and thermal route results

Canonical electrical series output is present for steady, 1-step, 7-step, and
HYPRE CPU/GPU cases. Existing electrical parity remains PASS, including the
nonzero 7-step pulse. HYPRE CPU/GPU scalar and vector parity also remains PASS.

Temperature continuity remains PASS: conformal shared-node jumps are zero and
Mortar nearest-coordinate jumps are at approximately the 1e-11 K scale.

The heat-flux checker was corrected to use the full tetrahedral 3D temperature
gradient and an outward normal derived from the parent tetrahedron centroid;
input triangle ordering is no longer trusted. A synthetic two-material
constant-flux test passes, including conductivity mismatch and reversed face
ordering. On the real 7-step result, all three interfaces have equal area and
facet counts, and left/right normals are consistently opposite:

| interface | normalized imbalance | status |
| --- | ---: | --- |
| Membrane/TES | 0.954546 | FAIL |
| TES/Stycast | 1.000992 | FAIL |
| Stycast/absorber | 0.157996, below informative flux scale | NOT_INFORMATIVE |

Therefore the normal-orientation and fixed-sign post-processing bugs are
**CLOSED**. Shared-face double counting is not evident. The remaining issue is
an open discrete elemental-flux reconstruction / physical-consistency gate;
it is not converted to PASS merely because route-to-route values are similar.

## Production readiness

| gate | status |
| --- | --- |
| electrical parity | CLOSED / PASS |
| HYPRE CPU/GPU smoke parity | CLOSED / PASS |
| production topology | CLOSED / PASS after regeneration |
| Stycast ideal/OCC/FEM convergence | CLOSED / PASS |
| real-case heat-flux validation | OPEN / FAIL |
| production HYPRE CPU/GPU benchmark | NOT READY / NOT RUN |
| production pulse | NOT READY / NOT RUN |

The GPU benchmark remains intentionally blocked until the heat-flux gate is
closed with a validated solver-native or converged flux diagnostic. No full
pulse was started.

## Artifacts

Key reports are in `artifacts/phase20_conformal/`:

- `heat_flux_diagnosis.json`
- `seven_step_physical_parity_flux_oriented.json`
- `production_topology_diagnosis.json`
- `production_current_conformal_interfaces.json`
- `production_geometry_volume.json`
- `stycast_convergence_report.json`
- `stycast_geometry_coarse.json`
- `stycast_geometry_medium.json`
- `stycast_geometry_fine.json`
- `electrical_parity.json`
- `production_benchmark.json`
- `production_pulse_waveform.json`

## Decision

**CONTINUE.** The production topology and geometry-understanding blockers are
closed. The conformal route is not yet ready for production GPU benchmarking
because the real-case heat-flux correctness gate remains open.
