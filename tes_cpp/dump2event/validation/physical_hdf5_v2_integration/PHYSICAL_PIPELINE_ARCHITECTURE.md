# Experimental physical HDF5 v2 pipeline

```text
raw Dumpall
  ├─ dump2event_physical_reconstruct (C++ parser + validated Pb ledger) ─┐
  └─ dumpall_transport_export (C++ parser; event and track CSV streams) ──┤
                                                                         └─ build_physical_hdf5_v2.py → EventV2Writer → atomic HDF5 v2
```

The programs read the Dumpall; they do not write it. Intermediate CSV rows are streamed event-by-event by Python and are not expanded for the entire run in memory. The measured 10,000-history conversion held one event's rows at a time, wrote 1,954,742 deposit records, and used about 301 MiB peak working set for the monitored Python process.

The only trusted reconstructors are the current Pb/PHITS 3.341 branches already in `PhysicalDepositLedger`: EGS5 continuous electron loss, boundary loss, electron cutoff, dead secondary electrons, and strict exact-zero rows. Unsupported rows become HDF5 `UNRESOLVED` with invalid energy. There is no legacy-formula fallback.

HDF5 `SEGMENT` start/end are raw transport endpoints. `score_xyz_cm` is left absent for segments because a representative score point is not independently derivable from these raw fields. Dead-electron `POINT` positions use raw secondary coordinates, which match the internal point oracle in the 100-history comparison. Zero rows retain a raw transport point but have zero energy. Parent-child links are unavailable and remain invalid. All event absorption labels remain `undetermined`.

This utility is explicit research tooling, not a CLI mode in normal `dump2event`, and does not add TES block allocation or Pulse conversion.
