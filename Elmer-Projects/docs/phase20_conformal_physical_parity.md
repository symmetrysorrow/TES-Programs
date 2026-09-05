# Phase20 conformal / Mortar physical parity

## Scope and worktree

This validation is on branch `phase20-conformal-physical-parity`, in the
separate worktree `D:/github/TES-Programs/Elmer-Projects-physical-parity`,
based on `origin/main` `8dd95cbd`. The earlier native-probe worktree was not
modified.

The existing Mortar + Block-Schur route remains the control. The new route is
conformal/shared-node, with the same physical mesh generated once and a
Mortar control mesh derived from it by duplicating only interface node IDs.

## Mesh and topology

The conformal parity mesh has 26,693 nodes and 132,357 tetrahedra. The control
mesh duplicates 79 Membrane/TES, 51 TES/Stycast, and 22 Stycast/absorber
interface node IDs. Coordinates, element geometry, body volumes, and contact
facets are otherwise preserved. The post-ElmerGrid conformal topology gate
passes on all three interfaces with zero-volume count zero.

The common mesh still under-resolves the literal 498 um Stycast cylinder:

- analytic volume: `3.8956377223e-12 m3`
- mesh volume: `3.7613195524e-12 m3`
- relative error: `-3.447912%`
- TES/Stycast contact area error: `-1.812701%`
- Stycast/absorber contact area error: `-6.451072%`

The convergence artifact records a regenerated base and nominal 2x level. They
are byte-identical in the resulting mesh despite different target values, so
they do not constitute a valid convergence series. Enabling contact-local
refinement at the requested scale caused the current Gmsh run to exit with a
native failure. This gate remains open.

## Electrical output and parity

The direct custom HeatSolve route is required for the fully coupled inner
circuit. The stock HYPRE CPU/GPU installation does not contain that custom
hook, so HYPRE cases use the portable external `circuit_parallel` UDF with one
external circuit update per solve. `run.py` now normalizes solver iteration CSV
or state output to the canonical schema:

`time_s,time_step,nonlinear_iter,tes_temperature_K,tes_current_A,tes_resistance_ohm,tes_power_W,bias_current_A,shunt_resistance_ohm`

Canonical output is now present for steady, 1-step, 7-step, and HYPRE CPU/GPU
cases. The comparison artifact reports PASS for all four comparisons:

- steady Mortar/conformal: max relative T/I/R/P =
  `3.72e-6 / 1.24e-4 / 1.69e-4 / 7.83e-5`
- 1-step Mortar/conformal: `3.20e-6 / 1.02e-4 / 2.42e-4 / 3.81e-5`
- 7-step Mortar/conformal: `8.70e-6 / 2.84e-4 / 4.37e-4 / 1.42e-4`
- HYPRE CPU/GPU at 1e-7: `1.83e-7 / 7.49e-6 / 9.37e-6 / 2.24e-7`

The nonzero 7-step pulse is electrically reproducible. Its Mortar-control
current amplitude is `4.08187e-6 A` (2.166% of baseline); the conformal
waveform differs by at most `5.47e-8 A`.

## Thermal continuity and heat flow

Temperature continuity passes: conformal shared-node jumps are exactly zero,
and Mortar nearest-coordinate jumps are at the 1e-11 K or lower scale in the
7-step result.

The current tetrahedral z-gradient heat-flux evaluator does not pass a
conservation gate on the nonzero pulse:

| interface | flux estimate | status |
| --- | --- | --- |
| Membrane/TES | `qL=-1.22618e-10 W`, `qR=2.69764e-9 W`, imbalance `0.954546` | FAIL |
| TES/Stycast | `qL=7.32324e-10 W`, `qR=7.26530e-13 W`, imbalance `1.000992` | FAIL |
| Stycast/absorber | magnitude below `1e-12 W` | NOT_INFORMATIVE |

These are reported as failures, not converted into a pass by route similarity.
The estimator may need a finite-element boundary-flux implementation, but
until that is independently validated the heat-flow gate remains blocking.

## HYPRE CPU/GPU

On the 26,693-node parity mesh, the current external-circuit HYPRE runs all
finish with three nonlinear iterations at 1e-7:

| run | HYPRE solve-time sum | total wall | result |
| --- | ---: | ---: | --- |
| CPU | `0.265716 s` | `4.19 s` | ALL DONE |
| GPU | `0.600553 s` | `4.36 s` | ALL DONE; GPU migration messages |

Vector parity is `max_abs=5.65e-6`, relative L2 `1.05e-6`. This mesh is too
small to claim GPU benefit; GPU is slightly slower. The current tolerance
sweep at 1e-6, 1e-7, and 1e-8 completes, but the 1e-8 result is only
comparable within the revised external-circuit route, not with the earlier
inactive-inner-circuit run.

## Production-size gates

An available 90,872-node / 459,683-tetra production candidate was checked
against the current branch. It fails the shared-node checker: the
Membrane/TES left surface is missing, TES/Stycast node IDs and coordinates do
not match, and Stycast/absorber partitioning differs. It is therefore not
used for benchmarking. No production-size CPU/GPU benchmark or production
pulse waveform is claimed; both are recorded as `NOT_RUN` with follow-up
requirements in the machine-readable artifacts.

## Artifacts and decision

Key reports are in `artifacts/phase20_conformal/`:

- `electrical_parity.json`
- `seven_step_physical_parity_nonzero.json`
- `pulse_waveform_parity.json`
- `stycast_mesh_convergence.json`
- `hypre_tolerance_study.json`
- `production_candidate_conformal_interfaces.json`
- `production_benchmark.json`
- `production_pulse_waveform.json`
- `mortar_control_provenance.json`

Focused tests cover the electrical-series normalizer and conformal/Mortar
parity. Unrelated legacy failures involving absent hybrid-prism and old
HeatSolve inputs were not expanded in this phase.

## Decision

**CONTINUE.** Do not PROMOTE. The conformal/shared-node implementation and
electrical/HYPRE route are technically viable, but production acceptance is
blocked by the heat-flux conservation failure, missing valid Stycast
convergence, and the invalid production candidate. The existing Mortar route
remains the fallback until those gates are closed.
