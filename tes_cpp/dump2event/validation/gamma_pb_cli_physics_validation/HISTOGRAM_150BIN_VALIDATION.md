# PHITS `[t-deposit]` 150-bin comparison

The run tally is `mesh=reg`, `reg=1000`, `part=all`, `deposit=0`, `unit=3`, `output=deposit`, `e-type=2`, 150 bins from 0 to 1.5 MeV (0.01 MeV wide). PHITS unit-3 counts are compared with per-history sums read back from HDF5, using the same run histories. Zero-deposit histories are excluded from the nonzero distribution; events above 1.5 MeV are reported as overflow instead of forced into the last bin.

| Measure | 0.511 MeV / 100 | 2 MeV / 100 | 0.511 MeV / 1,000 HDF5 | 2 MeV / 1,000 raw ledger |
|---|---:|---:|---:|---:|
| Zero-deposit histories | 2 | 9 | 27 | 87 |
| Underflow | 0 | 0 | 0 | 0 |
| Overflow >1.5 MeV | 0 | 19 | 0 | 275 |
| In-range PHITS / reconstruction counts | 98 / 98 | 72 / 72 | 973 / 973 | 638 / 638 |
| Absolute L1 bin-count difference | 0 | 0 | 0 | 0 |
| Maximum bin-count difference | 0 | 0 | 0 | 0 |
| Mean reconstructed event deposit (MeV, all histories) | 0.3794091954185283 | 0.9830887752521569 | 0.391834669051392 | 1.059594627322437 |
| Sample SD / standard error (MeV) | 0.18942508637608782 / 0.018942508637608782 | 0.6145946314727797 / 0.06145946314727797 | 0.17428499658286517 / 0.005511375511965169 | 0.6249594592643488 / 0.019762953365425606 |
| Peak bin | 0.51–0.52 MeV | 0.97–0.98 MeV | 0.51–0.52 MeV | 0.97–0.98 MeV |
| Maximum reconstructed event deposit (MeV) | 0.510995688 | 1.99997893 | 0.5109956880001679 | 1.9999808900003644 |

The 2 MeV / 1,000 comparison uses per-event sums from the reconstructed raw ledger, not HDF5: one below-guard secondary remains unresolved and the CLI withheld HDF5 publication. The ledger's supported rows exactly match the PHITS histogram counts, with 87 zero histories and 275 overflow histories. The bin-count equality is a same-run tally check. Event means include overflow; no comparison treats histogram midpoint counts as exact energy sums. These results do not imply complete absorption or pointwise scorer equivalence outside the compared region and settings.
