# HDF5 v2 implementation status

## Implemented in this change

- Added the experimental Python `EventV2Writer` and flat-offset v2 reader.
- Dataset groups are `/events`, `/particles`, `/tracks`, `/deposits`, and
  `/provenance`; fields and units are specified in `HDF5_V2_SCHEMA.md`.
- POINT, SEGMENT, and UNRESOLVED remain distinguishable. Segment endpoints
  and the PHITS representative score coordinate are separate.
- Unknown energy is returned as `None` with `energy_valid=false`; it is not
  scored as zero. Known negative, NaN, or infinite deposit values are rejected.
- Event, particle, and addition identifiers are checked for duplicates.
- The writer streams in batches of 256 events to extendible datasets, writes
  to a unique temporary file, validates offsets, and atomically publishes only
  after successful close. Failed conversion preserves a preexisting target.
- The reader dispatches by `format` and `format_version`, supports the Python
  flat v1 and legacy C++ hierarchy, and rejects unknown versions. v2
  `iter_events()` yields metadata/ranges; `iter_deposits()` yields typed rows.
- `Pulse_Ms()` preflights input before creating pulse output and rejects v2.
  Segment allocation has not been implemented in TES block space.

## Verification performed

- Focused Python tests: 6 passed after adding v2 and legacy compatibility
  coverage; comparator tests requiring unavailable external fixtures were
  skipped in the larger run.
- Relevant Python regression command:
  `PYTHONPATH=tes_cpp/python python -m pytest tes_cpp/tests tes_cpp/dump2event/tests -q`
  — 19 passed, 10 skipped.
- MSVC 19.51 CMake build: succeeded for the existing converter and physical
  ledger targets. CTest: 1/1 passed. The HDF5 schema code is Python and was
  covered by pytest.
- Synthetic write benchmark, 10,000 events / 20,000 tracks / 10,000 POINT
  deposits: 4.451698 s and 2,210,984 bytes on the current Windows workspace
  environment (about 2,246 events/s). This measures writer overhead and
  compression for the synthetic fixture, not a PHITS Dumpall conversion.
- Existing tracked PHITS validation reports remain: 38,077/38,077 additions
  in 100 histories and 201,971/201,971 additions in the independent 1,000
  history sample, maximum per-history difference
  `7.815970093361102e-14 MeV`, 150-bin histogram L1 difference 0. These were
  not rerun by the HDF5 storage change.

## Not implemented or not validated

- The default C++ converter still writes its legacy event representation.
  The validated `PhysicalDepositLedger` is still a separate CSV-producing
  research executable. It has not been connected to a complete streaming
  Dumpall-to-v2 pipeline; no opt-in physics CLI or HDF5 physical conversion
  was added here.
- Accordingly, the existing physical counts do not establish that all PHITS
  events or all photon/positron/weight cases are supported by v2.
- The v2 writer has not been driven from the real 10,000-history Dumpall or
  same-run PHITS internal addition IDs. No new PHITS comparison, deposit
  position match, or segment endpoint match was executed.
- HDF5 v2 is not ready for the everyday converter default or
  `Pulse_Ms()`. Keep the physical model opt-in and experimental until the
  source-ledger integration and same-addition position/energy checks exist.
- No change was made to the legacy `E_deposit` formula, the PR #6 source, the
  original PHITS input, the installed PHITS tree, or the source Dumpall.
