# Same-run internal addition reconciliation

Each 100-history PHITS diagnostic run emitted (1) a losslessly serialized raw Dumpall, (2) monotone raw `DUMP` record sequence numbers, (3) `TALLYSEQ` rows, and (4) paired PRE/POST internal tally traces. `compare_hdf5_to_phits_trace.py` joins by raw source record, exact sequence, and repeated-record occurrence; it does not infer an identity from histogram ordering.

| Quantity | 0.511 MeV, 100 | 2 MeV, 100 | 0.511 MeV, 1,000 | 2 MeV, 1,000 raw-only |
|---|---:|---:|---:|---:|
| PHITS PRE/POST additions | 250,938 | 549,639 | 2,541,564 | 5,978,912 |
| HDF5 deposits | 250,938 | 549,639 | 2,541,564 | not published |
| Raw ledger rows: reconstructed / unresolved | 250,938 / 0 | 549,639 / 0 | 2,541,564 / 0 | 5,978,912 / 1 |
| Raw transport records | 250,805 | 549,574 | 2,540,373 | 5,977,725 |
| Event / particle / NCOL / region mismatches | 0 | 0 | 0 | 0 |
| Terminal JCOLL carry-forward (NCOL 10/11) | not separated | not separated | not separated | 120,354 |
| Raw ledger vs PHITS `tlv`, max absolute error | — | — | — | 4.683146231920787e-14 MeV |
| Raw ledger vs PHITS `tlv`, sum absolute error | — | — | — | 7.540625087332713e-10 MeV |
| Rows with raw-ledger-vs-`tlv` error > 1e-12 MeV | — | — | — | 0 |
| Total reconstructed / PHITS internal energy | — | — | — | 1,059.594627322437 / 1,059.5946273224363 MeV |
| Difference of total energies | — | — | — | 6.82e-11 MeV |
| Event-summed maximum absolute difference | 1.1268763699945339e-13 MeV | 3.4794389591752406e-13 MeV | 1.1268763699945339e-13 MeV | 6.756817327868703e-13 MeV |
| Region / position mismatch | 0 / 0 | 0 / 0 | 0 / 0; max coordinate diff 0 cm | region 0; full position comparison unavailable without HDF5 |
| Negative / unresolved / unsupported | 0 / 0 / 0 | 0 / 0 / 0 | 0 / 0 / 0 | 0 / 1 / 1 |

The 100-history and 0.511 MeV/1,000 HDF5 joins use the same-run raw DUMP sequence sidecar and PHITS TALLYSEQ/PRE/POST traces. The 2 MeV/1,000 HDF5 was intentionally withheld because one candidate remains unresolved. Its 5,978,912 reconstructed raw-ledger rows were compared in order against the 5,978,912 internal additions; all event, particle, NCOL, and region identities matched. The 120,354 NCOL=10/11 JCOLL differences reflect PHITS carrying the originating JCOLL while the Dumpall terminal header stores zero. The energy and event-sum values above are from the reproducible raw-only comparison script, not an HDF5 result. Large raw/traces are intentionally not committed.
