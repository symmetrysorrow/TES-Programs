# Measured performance and output sizes

| Work | Histories | Wall time | Input/output size |
|---|---:|---:|---|
| PHITS diagnostic run R1000 | 1,000 | 76.74 s (PHITS-reported CPU time) | Dumpall 277,977,602 bytes; internal trace 164,150,585 bytes |
| Normal CLI HDF5 v2 R1000 | 1,000 | 28.72 s | HDF5 58,865,989 bytes |
| Normal CLI HDF5 v2 R100 | 100 | 3.65 s | HDF5 11,523,587 bytes |

These are measurements from this environment and these exact runs. Row-level diagnostic traces and raw files remain local. No 10,000-history HDF5 conversion was run in this turn.
