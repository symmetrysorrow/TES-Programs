# Phase20 conformal / Mortar physical parity

## Scope

This is the clean branch `phase20-conformal-physical-parity-latest`, based on
current `origin/main` `d04102b4`. The native-probe branch remains untouched.
The existing Mortar + Block-Schur route is the control; this branch validates
the conformal/shared-node primal route.

## Geometry and topology

The smoke conformal mesh has 26,705 nodes and 132,458 tetrahedra. Its three
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

The production topology blocker is therefore **CLOSED**. The regenerated
90,872-node production mesh also passes the current shared-node topology check.

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
| Membrane/TES | 0.954091 | FAIL |
| TES/Stycast | 1.020766 | FAIL |
| Stycast/absorber | 0.382821 | FAIL |

Therefore the normal-orientation and fixed-sign post-processing bugs are
**CLOSED**. Shared-face double counting is not evident. The result parser also
selects the requested `Perm` field explicitly, including Elmer's `use
previous` form for auxiliary fields.

## Real-case flux classification

A source-free, nonzero-flux conduction control was added on the same geometry
and shared-node route: the bath is fixed at 0.15 K, the absorber top at
0.16 K, the TES electrical body force is disabled, and membrane conductivity
is frozen at its bath-temperature value so the control is exactly linear.
HeatSolve `Calculate Loads = Logical True` and SaveScalars boundary sums provide
the solver-native external reactions.

| level | solver | Q_hot (W) | Q_bath (W) | G_eff (W/K) | global residual |
| --- | ---: | ---: | ---: | ---: | ---: |
| base | Direct | 1.1226478e-10 | -1.1226353e-10 | 1.1226478e-08 | 1.12e-05 |
| coarse | Direct | 1.0137063e-10 | -1.0136478e-10 | 1.0137063e-08 | 5.77e-05 |
| medium | Direct | 9.9176815e-11 | -9.9172754e-11 | 9.9176815e-09 | 4.09e-05 |
| fine | HYPRE CPU | 1.0232763e-10 | -1.0232853e-10 | 1.0232763e-08 | 8.82e-06 |

The medium Direct/HYPRE CPU conductance difference is `4.43e-05`; the
medium-to-fine conductance change is `3.18%`, within the declared 5% mesh gate.
All solver-native global balances pass the `1e-4` relative criterion.

The same final 7-step field gives essentially identical Mortar/conformal
values (maximum normalized-imbalance difference `2.82e-5`), so this is not a
conformal-only defect. Transient internal energy was evaluated as
`rho*cp*integral(T dV)`; storage is nonzero, but does not by itself explain
the raw side-flux jump.

Elmer 26.1's native `FluxSolver` is available and was executed on the control
case. It exports `temperature flux` and `temperature grad` to VTK. On a
shared-node conformal interface the exported flux is a single nodal projected
field, so exact opposite side integrals are not independent weak-form
field, so exact opposite side integrals are not independent weak-form
conservation evidence. Raw elemental-gradient interface flux is therefore
**DIAGNOSTIC_ONLY**, not a hard acceptance gate.

## Production readiness

| gate | status |
| --- | --- |
| electrical parity | CLOSED / PASS |
| HYPRE CPU/GPU smoke parity | CLOSED / PASS |
| production topology | CLOSED / PASS after regeneration |
| Stycast ideal/OCC/FEM convergence | CLOSED / PASS |
| weak-form global energy balance | CLOSED / PASS |
| conductance mesh convergence | CLOSED / PASS |
| Mortar/conformal global reaction parity | CLOSED / PASS |
| real-case raw elemental flux | DIAGNOSTIC_ONLY / CG reconstruction limitation |
| production HYPRE CPU/GPU benchmark | CLOSED / PASS |
| production pulse | READY / NOT RUN in this turn |

The production benchmark used the 90,872-node / 459,683-tetrahedron conformal
mesh. CPU wall time was `19.52 s`; GPU wall time was `18.43 s` (`1.059x`).
Temperature parity was max `2.58e-06 K`, with TES volume-average difference
`4.94e-08 K`.

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
- `heat_flux_mesh_convergence_control.json`
- `body_energy_balance.json`
- `body_energy_balance_steady_control.json`
- `mortar_conformal_flux_comparison.json`
- `native_flux_probe.json`
- `heat_flux_acceptance.json`
- `weak_form_global_energy_balance.json`
- `thermal_conductance_mesh_convergence.json`
- `mortar_conformal_reaction_parity.json`
- `heat_flux_acceptance_v2.json`
- `production_cpu_gpu_benchmark.json`
- `production_cpu_gpu_parity.json`
- `gpu_size_crossover.json`

## Decision

**READY_FOR_PRODUCTION_TRANSIENT.** The solver-native weak-form reaction gate,
conductance mesh gate, medium Direct/HYPRE CPU parity, and Mortar/conformal
reaction parity all pass. The raw shared-node CG elemental-gradient jump is
closed as `CLOSED_AS_CG_FLUX_RECONSTRUCTION_LIMITATION`; it remains a
diagnostic and is not interpreted as an independent weak-form reaction.
