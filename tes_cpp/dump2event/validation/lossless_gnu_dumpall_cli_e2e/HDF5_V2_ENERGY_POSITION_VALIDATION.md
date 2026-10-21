# HDF5 v2 measured energy, position, and count audit

R100 and R1000 HDF5 files were produced through the normal CLI and reread using the repository's HDF5 v2 reader. Same-run row-level comparison used the scripts in `scripts/` and PHITS PRE/POST plus DUMP/TALLYSEQ sidecars.

| Metric | R100 | R1000 |
|---|---:|---:|
| Events | 100 | 1,000 |
| HDF5 deposit additions | 38,077 | 195,185 |
| `POINT` / `SEGMENT` | 160 / 37,917 | 1,267 / 193,918 |
| Exact-zero addition rows | 130 | 1,141 |
| Negative / unresolved rows | 0 / 0 | 0 / 0 |
| PHITS internal total | 6.929495555108982 / 34.510355862665456 MeV | |
| HDF5 total by event group | 6.929495555108977 / 34.510355862665435 MeV | |
| Per-addition max absolute error vs internal increment | 5.80937208061183e-15 / 5.818913059729702e-15 MeV | |
| Per-event sum max absolute error | 3.26405569239796e-14 / 5.129230373768223e-14 MeV | |
| Region mismatches | 0 | 0 |
| Position rows matched / max coordinate difference | 38,077 / 0 cm | 195,185 / 0 cm |

The 1000-history event distribution has mean `0.03451035586266543 MeV`, sample SD `0.17010301097575625 MeV`, standard error `0.005379129515360106 MeV`, and maximum `1.3319859860001428 MeV`. There are 944 zero-total histories and 56 positive-deposit histories. The 100-history run has 91 zero-total histories and 9 positive-deposit histories.

All HDF5 event records are marked complete for this candidate scoring slice. The HDF5 provenance still says `transport_physics_completeness=undetermined`, `all_absorption=undetermined`, and `energy_balance_status=not assessed`; “complete” here only means no unsupported rows were emitted under this experimental algorithm. No claim is made about photon/positron/general-material coverage.
