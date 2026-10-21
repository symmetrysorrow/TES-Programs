# Validation provenance audit

The evidence contains three distinct run artifacts. They must not be silently treated as one run.

| Artifact | Run identity | Dumpall size / SHA-256 | Input / executable evidence | Relation |
|---|---|---|---|---|
| Earlier 10k source dataset | `overnight_physics_reconstruction/runs/n10000_independent_seed` | 2,783,824,000 bytes / `d254b3afd13230be84ad8f733c8aaa89811724f832a506a0d5b421d051c06b1c` | Input has `maxcas=10000`, `maxbch=1`, `negs=1`, `e-mode=2`, cutoff 0.001 MeV, and seed `010110…1011`; PHITS 3.341. | This is the raw file from which the earlier HDF5 candidate was produced. |
| 2026-10-21 baseline rerun | `D:\temp\positron_10000_closure_20261021\baseline` | 2,736,869,598 bytes / `c11c51299c6c77a7e8a531a640ed72b31c8455fb2d0a17f2b03d3716a2815060` | `phits_baseline.exe`, SHA-256 `80496c548131bf214601db09005d305e06c83afc116dDA92a9b65790d700a5c2` (case-insensitive); same 64-bit input seed and physical settings. | Its raw Dumpall is byte-identical to the diagnostic rerun below, not to the earlier raw file. |
| 2026-10-21 full diagnostic rerun | `D:\temp\positron_10000_closure_20261021\diagnostic_full` | 2,736,869,598 bytes / same `c11c51299…` | `phits_WinGfort.exe`, SHA-256 `0477c7c6a0849ad1ac6e5ac97b149f55bdaf623c4ddc3eff407b738bead8c51b`; input template SHA-256 `7401c49903c41317825af7eed2a8c239033ea7505514cecebe0f1524459f007e`. | Baseline and diagnostic `[t-deposit]` files are byte-identical, SHA-256 `0d87d8b7041c54cfd10e610ac0653ff1eda1844a31ca0691126ce82ff0c5f110`; the internal trace has 1,954,742 PRE rows. |

The two input files have the same transport settings and seed; their textual differences are the `file(15)` path and trailing newline. Their PHITS executables differ in hash/build. Thus the Dumpall hash discrepancy between the earlier run and the October 21 reruns remains **unresolved**; identical seeds and tally files do not prove raw-record identity. The October 21 baseline and diagnostic Dumpall hashes match each other, which supports diagnostic non-perturbation for that pair.

The structured MSVC parser accepts the earlier `d254…` raw Dumpall: 1,983,461 records, 10,000 histories, terminal NCOL=2. It rejects the October 21 GNU-built diagnostic Dumpall at physical line 5: an NCOL=4 particle tuple has 11 combined numeric values where 8 are required, consistent with the known 80-column `cpcasdump` truncation. That file is not repaired by guessing. Therefore the prior HDF5/internal match is not accepted as a same-run raw/internal proof until the exact files used for the join are re-audited or the new diagnostic raw is regenerated losslessly.

Run path, file size, mtime, hashes, and exact local commands are recorded in the ignored local audit file `private/PROVENANCE_LOCAL.csv` when available. Private PHITS source, internal traces, raw Dumpall, HATCH table, and HDF5 outputs are not included in this public report or Git.
