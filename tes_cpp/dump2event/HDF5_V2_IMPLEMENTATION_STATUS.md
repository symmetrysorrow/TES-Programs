# HDF5 v2 implementation status

Status: experimental storage API with a real-Dumpall diagnostic pipeline. The
ordinary C++ converter still uses its legacy output by default; no production
Pb mode has been added to that CLI.

## Implemented

- `EventV2Writer` writes the flat-offset v2 schema in batches to a private
  temporary file and atomically publishes only after successful close.
- `append_event()` validates and stages the full event before mutating buffers,
  offsets, or global identifiers. Caught validation exceptions leave no event
  rows behind. A storage failure poisons the temporary writer and blocks later
  appends/publication.
- The reader supports v2, Python v1, and the legacy C++ hierarchy. The v2 API
  preserves POINT, SEGMENT, and UNRESOLVED; unknown energy is never read as a
  known zero.
- HDF5 provenance records the actual Dumpall, PHITS input, stopping table, and
  canonical configuration SHA-256 values for the new experiment.
- An opt-in validation utility streams raw PHITS records through the existing
  Pb physical ledger and raw transport exporter into HDF5 v2. It is separate
  from the ordinary converter and reports unsupported rows as UNRESOLVED.
- `Pulse_Ms()` still preflights and rejects v2; TES block allocation and pulse
  generation remain outside this format conversion.

## Real PHITS checks

- 100 histories: 38,077/38,077 HDF5 readback deposits match PHITS internal
  addition IDs, history, particle, region, energy (maximum error
  `5.7976626971489864e-15 MeV`), and checked positive position. Missing/extra
  additions and checked position mismatches: 0.
- Independent 1,000 histories: 201,971 HDF5 deposits; 150-bin `[t-deposit]`
  histogram L1 difference 0. This is a histogram-level check, not per-addition
  internal matching.
- 10,000 histories: 1,954,742 ledger rows, of which 342 are unresolved. The
  HDF5 readback histogram differs in 13 bins (L1=18). This is partial and
  does not pass physical completeness.

## Limits

- The normal `dump2event` CLI does not yet expose `--physics pb-egs5` or write
  HDF5 v2. Parent-child linkage and input-particle metadata are not reconstructed.
- `complete` describes whether the configured ledger slice has unresolved
  rows; it does not prove complete transport or energy conservation.
  `all_absorption` remains `undetermined` for every written event.
- Continuous loss endpoints are preserved as SEGMENT endpoints. Their PHITS
  representative score coordinate is not inferred from them and remains
  absent. The 100-history test checks interval endpoints and validated point
  positions separately.
- Variable weights, positron annihilation, other photon branches, and the 342
  unsupported 10,000-history rows remain unresolved. No negative clipping,
  normalization, or legacy-formula fallback is used.
- HDF5 files, Dumpalls, PHITS source, stopping tables, and internal trace rows
  remain local and are not committed.

See [`validation/physical_hdf5_v2_integration/`](validation/physical_hdf5_v2_integration/)
for the exact counts and run comparisons.

## 2026-10-21 same-run positron closure update

An isolated PHITS 3.341 diagnostic build and the same 10,000-history Pb run
closed the earlier 342 unresolved additions and identified 1,460 positron
continuous-loss additions that had been misclassified as exact zero. The
separate experimental ledger now uses PHITS EGS5 positron PDedx coefficients
for primary transport and QCLUSTS(7) kinetic residual for the validated dead
positron secondary branch. A newly written HDF5 v2 has 1,954,742/1,954,742
known nonnegative deposits. Same-run internal addition comparison has zero
unresolved, energy, region, or checked positive-position mismatches; maximum
addition error is `5.865533753146579e-15 MeV`. The 150-bin PHITS history
histogram has L1 count difference 0 (the PHITS printed decimal residual is
`1.554312234475219e-14`).

This result is scoped to that PHITS 3.341 seed, Pb region 1000, medium 2,
unit weights, explicit 1 keV cutoffs and supplied matching stopping table.
Variable weights, other materials, photon-process-specific raw-only ledgers,
and independent-seed internal-addition validation remain unverified. This
does not change the default converter or PR #6.
