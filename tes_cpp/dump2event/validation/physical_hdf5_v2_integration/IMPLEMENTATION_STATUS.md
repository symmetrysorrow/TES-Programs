# Physical Dumpall to HDF5 v2 integration status

Run date: 2026-10-21. Worktree base: `806e681f73e845b0e71d4bcf64a951fbf2356b8a`.

## Completed in this experiment

- `EventV2Writer.append_event()` now validates all event, particle, track, and deposit fields in local staging buffers before modifying writer buffers or ID sets. Buffer extension failures roll back the staged in-memory rows. An HDF5 storage failure poisons the private writer and prevents publication.
- Added transactional rejection/recovery tests and wired v1/v2 Python tests plus existing C++ tests into GitHub Actions.
- Added an explicit research-only raw-record exporter and a streaming Pb-ledger-to-HDF5 v2 utility. They are separate from the ordinary `dump2event` CLI and cannot silently alter its default behavior.
- Generated fresh HDF5 v2 outputs from the repaired 100-history, independent 1,000-history, and 10,000-history PHITS 3.341 Dumpalls. Outputs are in the ignored local `build/physdeposit-vs/` directory, not source control.
- HDF5 provenance stores exact Dumpall, PHITS input, stopping-table, and canonical physics-configuration SHA-256 values.
- Re-read each HDF5 v2 output. The 100-history file matches all 38,077 PHITS internal addition rows. The 1,000-history HDF5 reproduces the 150-bin PHITS `[t-deposit]` distribution exactly. The 10,000-history HDF5 contains 342 explicit unresolved rows and its histogram differs from PHITS; it is a partial diagnostic output.

## Not complete

- The research utility is not an option in the ordinary `dump2event` CLI. It invokes the existing dedicated experimental reconstruction executable.
- The track stream stores raw records and positions, but parent-child links and source-particle metadata are not reconstructed. `all_absorption` is always `undetermined` because complete energy balance and all transport branches have not been established.
- For 10,000 histories, 333 unsupported positron-parent scoring rows and 9 unsupported secondary rows remain unresolved. No score is substituted from the old parent-minus-secondary formula.
- There is no same-addition PHITS internal ledger available for the independent 1,000- or 10,000-history runs in the local validation inputs. Their comparisons are histogram-level only.
- No output HDF5, Dumpall, stopping table, PHITS source, or internal trace is committed.

See the companion reports for exact counts, hashes, error limits, and test results. The experimental scope is not approved for production or PR #6.
