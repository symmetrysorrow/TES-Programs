# Pb HATCH stopping table provenance

Historical manifest `D:\temp\physdeposit_publish_20261020\ledger.csv.manifest.json` records the expected private table path `D:\github\TES-Programs\dump2event_physics_validation\egs5_stopping_power_reconstruction\runs\n100_stopping_probe_v2\egs5_stopping_power_coefficients.csv` and SHA-256 `20fd86eb0b57483bd19f20a8388d2c1ecb0431a02fa99f0933e2c3c4748514e9`. The table file was not present at that path in the current workspace. Searches by exact filename under `D:\temp`, the research worktree, and the user's Documents/Downloads found no candidate. No candidate file was accepted based only on a similar name or table shape.

The old manifest identifies the source as a local PHITS 3.341 HATCH EGS5 Pb medium-2 probe and records a unit-weight, 0.001 MeV cutoff configuration. The coefficient file itself, row count, EDedx/PDedx values, HATCH setup transcript, and underlying PHITS executable hash are not available for re-verification here. No replacement table was generated. Stopping-table hash for this task: **not recovered**.

**Status:** required private input unavailable; physical CLI run blocked. Do not substitute a generic stopping-power table.
