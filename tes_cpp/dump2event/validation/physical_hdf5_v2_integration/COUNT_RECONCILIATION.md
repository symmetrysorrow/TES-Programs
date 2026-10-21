# Count reconciliation

| Dataset | Raw records (terminal marker excluded) | Histories | PHITS additions / reference level | HDF deposits | HDF zero / point / segment / unresolved |
|---|---:|---:|---|---:|---|
| 100 | 38,345 | 100 | 38,077 internal addition rows | 38,077 | 130 / 160 / 37,917 / 0 |
| 1,000 | 204,835 | 1,000 | histogram only; 58 positive histories | 201,971 | 1,163 / 1,306 / 200,665 / 0 |
| 10,000 | 1,983,461 | 10,000 | histogram only; 532 positive histories | 1,954,742 | 12,910 / 14,220 / 1,940,180 / 342 |

The 10,000-history table is the `1,983,461`-record file. The separate `2,036,435`-record input was not present in the inspected validation run directories and is not counted or inferred here. Raw record counts are parser records, not tally additions.

The 100 run splits into 37,190 continuous electron additions, 727 electron cutoff additions, 30 dead-electron secondary point additions, and 130 exact-zero additions. The 1,000 and 10,000 dump ledgers additionally include 8 and 83 EGS5 boundary additions respectively.
