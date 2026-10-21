# Real Dumpall end-to-end

## Existing raw audit

The earlier source dataset at `dump2event_physics_validation/overnight_physics_reconstruction/runs/n10000_independent_seed/dumpall.dat` was read by the structured parser: 1,983,461 records, 10,000 histories, terminal NCOL=2. Its public hash is recorded in `VALIDATION_PROVENANCE_AUDIT.md`.

The October 21 GNU diagnostic Dumpall has the same-size baseline and diagnostic files and matching SHA-256, but the current MSVC structured parser fails on the first NCOL=4 logical particle tuple (11 numeric values combined where schema requires 8). Source audit attributes the serialization loss to the GNU ASCII `cpcasdump` 80-character merge buffer. It is not repaired by inference.

The new ordinary CLI was exercised on this diagnostic file with the explicit Pb options. It exited nonzero when the transport exporter rejected physical line 5; the requested output HDF5 did not exist afterward. This confirms fail-closed dispatch for this malformed input, not successful conversion.

## Integrated CLI

A successful `dump2event --physics pb-egs5` conversion of 100, 1,000, or 10,000 real histories was **not run** in this task. The local private Pb HATCH coefficient CSV required by the tested engine was absent. Consequently no new integrated CLI HDF5 data, runtime, capacity, or memory measurement is reported. The historical HDF5 is not relabeled as an integrated CLI result.
