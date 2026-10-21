# PHITS `[t-deposit]` histogram reconciliation

Comparison script: `scripts/compare_tdeposit_histogram.py`. It parses the printed histogram edges and compares HDF5 event totals on the same edge set. It explicitly excludes exact-zero event totals because the positive-score PHITS output has no zero-score bin entries. `unit=3` values are per source; multiplying by the run's history count recovers integer history counts. Bins are `[e_i,e_(i+1))` with 0.01 MeV width, edges 0 to 1.5 MeV.

| Run | Histories | PHITS positive-history count | HDF5 positive-history count | L1 count difference | Max per-bin difference | Under/overflow |
|---|---:|---:|---:|---:|---:|---:|
| R100 | 100 | 9 | 9 | 0 | 0 | 0 / 0 |
| R1000 | 1,000 | 56 | 56 | 0 | 0 | 0 / 0 |

R1000 peak bin is 1.17–1.18 MeV. R100 peak bin is 1.33–1.34 MeV. The complete per-bin aggregate CSVs are local private artifacts and are not included in Git. This histogram result is an aggregate check; the stronger evidence is the same-run addition-level join in `RAW_INTERNAL_ADDITION_MAPPING.md`.
