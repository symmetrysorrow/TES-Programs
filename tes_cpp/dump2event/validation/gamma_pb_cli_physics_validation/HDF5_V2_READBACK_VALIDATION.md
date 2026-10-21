# HDF5 v2 CLI and readback validation

Both 100-history gamma datasets were converted with the normal `dump2event.exe` CLI using explicit `--physics pb-egs5 --hdf5-version 2 --physics-config ... --stopping-table ...`; `--normalize-deposits` was not used. The CLI published a v2 file only after all ledger rows were classified as reconstructed.

| Run | HDF5 events | particles | tracks | deposits | unresolved | negative | HDF5 bytes | Provenance checks |
|---|---:|---:|---:|---:|---:|---:|---:|---|
| 0.511 MeV / 100 | 100 | 4,955 | 250,805 | 250,938 | 0 | 0 | 75,482,132 | revision 2; CLI/HDF5 hashes equal |
| 2 MeV / 100 | 100 | 11,393 | 549,574 | 549,639 | 0 | 0 | 164,283,437 | revision 2; CLI/HDF5 hashes equal |
| 0.511 MeV / 1,000 | 1,000 | 50,406 | 2,540,373 | 2,541,564 | 0 | 0 | 763,709,732 | revision 2; raw/tally join complete, zero mismatches |
| 2 MeV / 1,000 | no published HDF5 | 5,977,725 raw tracks exported | — | — | 1 | 0 | no HDF5 | CLI exit 1; one sub-guard candidate retained unresolved |

The 2 MeV HDF5 exact addition join reports 549,639 matched rows. It has zero event/particle/region mismatches and no unmapped/overmatched raw sequence occurrences. The HDF5 files retain `complete_semantics` as candidate scoring completeness, `all_absorption=undetermined`, and `transport_physics_completeness=undetermined`; they do not claim all energy is absorbed in Pb.

CLI stdout JSON and `/provenance` attributes were compared with `scripts/verify_hdf5_provenance.py`. The compatibility alias equals the canonical resolved hash. Contract JSON, canonical settings, PHITS input, PHITS executable, converter, reconstruction helper, exporter, driver, Dumpall, and stopping-table digests have distinct meanings. The 1,000-history HDF5 v2 readback returned `hashes_match=true` with 1,000 events and 2,541,564 deposits.

The 2 MeV / 1,000 run's raw ledger still has an addition-level comparison: all 5,978,912 reconstructed rows match the PHITS internal increments within 4.69e-14 MeV, and the raw-ledger event histogram matches the PHITS unit=3 histogram. The one unresolved below-guard secondary means there is no complete HDF5 file and no HDF5 readback claim for that run.
