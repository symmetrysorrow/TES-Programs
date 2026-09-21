# No-mortar interface and COMSOL convergence study

## Main result

The same conformal mesh with mortar enabled changed the steady current by only -0.0153 uA.  The large mismatch was therefore not caused by the mortar flag alone.  Refining the no-mortar mesh moved the steady current from 191.673 uA (coarse) to 154.859 uA (fine), and to 144.764 uA at iteration 25 of the 20/40 um probe; COMSOL is 143.055 uA.  Spatial resolution in the TES/Stycast stack is the dominant demonstrated effect.

## Mesh and steady convergence

| case | nodes | tetrahedra | TES tets | Stycast tets | TES volume [m3] | final current [uA] | vs COMSOL |
|---|---:|---:|---:|---:|---:|---:|---:|
| coarse no-mortar | 26,705 | 132,458 | 508 | 317 | 4.000000e-14 | 191.673156 (+33.986%) |
| fine no-mortar | 90,872 | 459,683 | 1,615 | 901 | 4.000000e-14 | 154.858777 (+8.251%) |
| 20/40 um no-mortar probe | 538,253 | 2,826,834 | 6,494 | 3,831 | 4.000000e-14 | 144.710781 (+1.157%) |
| mortar reference (different mesh) | 63,172 | 312,997 | 7,929 | 2,865 | 4.000000e-14 | n/a |

The refine20 steady solve completed with exit 0 at nonlinear iteration 21 under the 1e-7 convergence setting and produced a restart `.result`. The earlier 25-iteration interrupted probe is retained separately as additional convergence evidence.

## Interface checks

- No-mortar temperature is topologically continuous because the interface uses shared nodes and therefore a single temperature degree of freedom. A numerical two-sided heat-flux jump was not measured because the stored runs have no VTU field output.
- Electric potential and current-density continuity are not applicable to this Elmer model: electrical behavior is represented by the lumped TES circuit, not a distributed electric PDE.
- The common COMSOL/Elmer circuit, material, pulse, and initial-state constants are already parity-checked in `nomortar_parity_diagnostic/diagnostic.json`; the mesh-dependent pulse discrete norm remains a separate discretization quantity.

## Transient convergence

The existing full waveform comparison remains confounded by different meshes: no-mortar fine has baseline 154.851 uA, peak drop 6.266 uA, and delay 0.100 ms; the mortar reference has 148.130 uA, 7.725 uA, and 0.500 ms; COMSOL has 143.055 uA, 7.775 uA, and 0.428 ms.

The completed refine20 short probe gives 144.711034 uA at 1 ns and 144.725127 uA at 101 ns after the event (about +0.000253 and +0.014346 uA relative to its steady 144.710781 uA). This confirms a stable, very small immediate response, but it does not determine the full 100 us-scale amplitude or time constant. A full transient refine20 waveform remains computationally expensive with MUMPS.

## Files

- `report.json` contains machine-readable mesh, steady, parity, and field-check status.
- `results/case_tes_steady_singlepixel_conformal_gpu_refine20/solver.log` and `.../case_tes_steady_singlepixel_conformal_gpu_refine20_iterations.csv` contain the probe evidence.
