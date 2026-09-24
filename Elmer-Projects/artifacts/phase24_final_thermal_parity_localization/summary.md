# Phase24 final thermal parity localization

Diagnostic-only campaign. Physics parameters, materials, TES law, circuit constants, bath BC, and production mesh were unchanged. Mesh generation was serialized and the new intermediate used a unique work directory.

## Fixed-power convergence

| case | G_eff (W/K) | historical ratio | interpretation |
|---|---:|---:|---|
| Historical | 1.878294461631e-08 | 1.000000 | reference |
| Refined mortar parent | 2.029899315000e-08 | 1.080714 | audited parent |
| Membrane/substrate intermediate (h=12.5 µm) | nan | nan | CPU/MUMPS if gate passed |
| Membrane/substrate control (h=10 µm) | 1.959713948000e-08 | 1.043348 | audited accepted control |

The attempted intermediate failed the pre-solver mesh gate and is excluded from the valid convergence series; no solver result is accepted for it.

An earlier h=20 µm exploratory attempt was also rejected after audit (Stycast upper trace 4579→1165 and bath faces 18592→18620); its 2.033979499e-8 W/K result is retained only in `rejected_h20_exploratory.json`.

## TES/Membrane gate

The previous TES/Membrane candidate is `invalid controlled mesh`. Its non-target statistics change (including the Stycast upper trace), so no TES/Membrane fixed-power result is accepted and its 2.6857e-8 W/K value remains confounded. `f_TM` is therefore not identifiable.

## Branch and trace diagnostics

`parallel_branch_heatflow.csv` records Q_i/G_i as unavailable: the existing native captures have no separate flux mask or mortar reaction for the conformal downstream Si1/SiNx and Si2/SiO2_2 paths. Splitting total P by area or assuming a series chain would be invalid. Trace slopes and differential segment resistances are materialized in `trace_temperature_slopes.csv` and `differential_resistance_breakdown.csv`.

The accepted Membrane/substrate control explains 46.29% of the parent-to-historical G_eff gap; 53.71% remains unexplained by valid one-factor tests. Existing resistance decomposition assigns 45.8% to TES/Membrane effective discrete response and 54.2% to the downstream network, but only the downstream share has a clean controlled mesh test.

## GO/NO-GO

Best valid thermal parity is `1.043348` (4.3348% high), outside the 2–3% GO band; full nonlinear CPU/MUMPS validation was skipped. Final current parity is therefore not measured. HYPRE/GPU: **NO-GO**.

Full ElmerSolver status: the audited parent and accepted h=10 µm control each have completed 0.95/1.00/1.05P CPU-native MUMPS captures with exit 0; the h=12.5 µm candidate was not launched because the pre-solver gate failed. Physics/material/TES law/circuit/bath and production mesh: unchanged.

See `mesh_statistics.csv`, `mesh_control_gate.csv`, `remaining_gap_accounting.json`, and the trace/branch CSVs for the machine-readable audit.
