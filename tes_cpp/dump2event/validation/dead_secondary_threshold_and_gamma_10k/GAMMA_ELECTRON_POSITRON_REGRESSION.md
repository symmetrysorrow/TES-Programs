# 511 keV / 2 MeV gamma回帰実測

ここでは既存のPHITS 3.341回帰実験の数値を固定する。これらは10,000履歴の独立seed検証の代替ではない。

| source / histories | Internal additions | HDF5 v2 rows | ID/event/particle/region/position mismatch | Negative / unresolved | 150-bin tally L1 |
|---|---:|---:|---:|---:|---:|
| 511 keV / 100 | 250,938 | 250,938 | 0 | 0 / 0 | 0 |
| 511 keV / 1,000 | 2,541,564 | 2,541,564 | 0 | 0 / 0 | 0 |
| 2 MeV / 100 | 549,639 | 549,639 | 0 | 0 | 0 |
| 2 MeV / 1,000 (new HDF5 readback) | 5,978,912 | 5,978,912 | 0 | 0 / 0 | 0 |

2 MeV / 1000 exact join: max addition energy error `4.680023729664029e-14 MeV`, total absolute error `5.133758390738264e-10 MeV`, max history sum error `6.756817327868703e-13 MeV`, and maximum score-coordinate error `0 cm`. HDF5 retains 4,797 exact-zero scored contributions and one independently proven `NOT_SCORED_BY_PHITS` candidate outside the deposit table. No unsupported or unresolved contributions remained for that run.

The 511 keV and earlier 2 MeV values are existing independent regression datasets described in prior validation reports. The new 10,000-history run below remains untested until PHITS, normal CLI, and same-run internal-trace joins all finish.
