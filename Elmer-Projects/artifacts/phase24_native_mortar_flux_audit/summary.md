# Phase24 native mortar heat-flow audit

## Scope

This is a diagnostic audit of the existing CPU MUMPS native full saddle-point captures. Production physics, materials, TES law, circuit constants, geometry, mesh topology, mortar formulation, timestep, and HYPRE/GPU settings were not changed.

## Native reaction definition

- The solver-native capture supplies the assembled constraint block `B`, the transpose block `C`, and the solved multiplier vector `lambda`.
- Native capture location: `../tools/elmer-hypre/src/fem/src/SolverUtils.F90`, `SolveWithLinearRestriction` / `Phase24CaptureFullSystem`; audit conversion location: `scripts/support/run_phase24_native_mortar_flux_audit.py`.
- The integrated reaction is `C lambda` summed over the primal DOFs belonging to the interface side; this is an actual constraint reaction, not a one-sided gradient reconstruction.
- Galerkin thermal projector weights carry m², multiplier values carry W/m², and the integrated primal reaction carries W.
- `B=slave-master`; slave and master reactions should cancel. The row-level records preserve row ID, DOFs, multiplier, both side reactions, sign, and units.

## Results

- Historical TES→membrane reaction: **3.203008493e-10 W**; Phase24 TES→membrane mortar rows: **0**.
- Phase24 mortar reaction is instead limited to TES↔Stycast and Stycast↔substrate rows; the remaining native bath heat is classified as conformal/shared-node or otherwise non-mortar path: **4.233820811e-10 W**.
- Historical native balance: Joule=3.203004762e-10 W, bath=3.203008670e-10 W, error=-3.908273862e-16 W (-1.220e-06).
- Phase24 native balance: Joule=4.232286181e-10 W, bath=4.233820809e-10 W, error=-1.534627568e-13 W (-3.626e-04).
- Matched-T (168.563176 mK) secant Q estimate: historical=3.203010517e-10 W, Phase24=4.746491027e-10 W; Phase24/historical=1.481884.
- Matched-power (3.203004762e-10 W) secant T estimate: historical=168.563142 mK, Phase24=162.526715 mK.
- A symmetric ±δT local derivative was not run; the recorded G_eff values are native-Q secants and are explicitly not a local derivative proof.

## Required conclusions

1. Native mortar reaction physical integration: **YES**, from native `C lambda`; row-level and aggregate artifacts are written.
2. Reaction-inclusive energy closure: **YES for global native Joule/bath closure**; the interface partition closes only for the mortar portion, with conformal/shared-node heat reported separately.
3. Same-T conductance: the native operating-point secant estimates Phase24 higher than historical; exact matched-T solve remains pending.
4. Same-P temperature: the secant estimate gives a lower Phase24 TES temperature, indicating higher effective transport in the captured native bath path.
5. Dominant path difference: historical TES→membrane is an explicit mortar path; Phase24 has no TES→membrane mortar rows and carries heat through a non-mortar/conformal path in this capture.
6. Area alone: do not claim area sufficiency; the interface-area CSV is provided, but the zero-vs-nonzero constraint topology is the larger structural difference.
7. Mortar topology: **YES, materially different**; per-interface counts and row norms are recorded.
8. 143→218 µA: **not quantitatively proven by this audit alone**. It establishes a large topology/path candidate, but the matched-point rows are secant estimates rather than new frozen-source solves.
9. Minimum next controlled fix: make the TES–membrane interface coupling identical between meshes (preserve the Phase24 physical mesh, add only the missing controlled interface representation), then rerun native MUMPS with ±δT and fixed-P diagnostics.
10. HYPRE/GPU: **NO-GO** until the controlled thermal-path comparison is completed.

## Artifacts

- `native_mortar_reaction_integrals.csv` — integrated reaction by case/interface.
- `native_mortar_reaction_rows.csv` — constraint row ID, primal DOFs/nodes, multiplier, side reactions, units, and sign.
- `native_heatflow_partition.csv` — reaction-inclusive path partition and native global closure.
- `matched_temperature_comparison.csv`, `matched_power_comparison.csv`, `effective_conductance_comparison.csv` — explicit secant diagnostic status.
- `interface_area_comparison.csv`, `mortar_topology_comparison.csv`, `energy_balance_closure.json` — geometry/topology and provenance.
