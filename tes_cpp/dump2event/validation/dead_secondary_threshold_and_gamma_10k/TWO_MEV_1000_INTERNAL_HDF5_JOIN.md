# 2 MeV / 1000履歴：PHITS内部加算とHDF5 v2の全件join

判定: **成功**。2026-10-21に実行したPHITS 3.341の同一runを用い、通常 `dump2event` CLIが書いたHDF5 v2を、実測された `internal_additions.csv` / TALLYSEQへ全件joinした。これはHDF5出力のreadback結果であり、raw-onlyの既存集計を代用していない。

## 来歴

- PHITS run: `D:\temp\gamma_pb_cli_physics_20261021\seq_n1000_2mev\`
- 元Dumpall: `gamma_2mev_thickPb_dumpall.dat`、8,457,662,268 bytes、SHA-256 `f84d0dd4624549bfcc1b54e705109aecc93c3b9a7ee345089d16d49baabc59dd`
- CLI HDF5: `H:\TES-Programs-validation\dead_secondary_threshold_and_gamma_10k\two_mev_1000\gamma_2mev_1000.h5`、SHA-256 `e90f120d548d3b128f6cac4b9068e8b749ffd013eb4997f51bd068bcf576d358`
- PHITS内部追加CSV: `internal_additions.csv`、5,978,912行
- Join結果: `internal_hdf5_join.json`
- Audit: `hdf5_counts.json`, `hdf5_tally_comparison.json`, `hdf5_tally_comparison.csv`

実行コマンドは `dump2event.exe <dumpall> <h5> --physics pb-egs5 --hdf5-version 2 --physics-config <pb_config.json> --stopping-table <egs5_stopping_power_coefficients.csv>`。標準CLI runのログは `cli.log`。PHITS tallyは `gamma_2mev_thickPb_tdeposit.out`、region 1000、unit=3、0–1.5 MeV、150 bins。

## 加算単位の照合

| 項目 | 結果 |
|---|---:|
| Internal additions / HDF5 deposit rows | 5,978,912 / 5,978,912 |
| Addition ID range / duplicate / missing | 1–5,978,912 / 0 / 0 |
| event, particle, region mismatch | 0, 0, 0 |
| position compared / unavailable / mismatch | 5,978,912 / 0 / 0 |
| Maximum position error | 0 cm |
| HDF5 energy vs TLV max abs error | 4.680023729664029e-14 MeV |
| Sum of absolute energy errors | 5.133758390738264e-10 MeV |
| Rows with error > 1e-12 MeV | 0 |
| PRE increment vs PHITS TLV maximum error | 1.1102230246251565e-16 MeV |
| HDF5 total / internal TLV total | 1059.594627322437 / 1059.594627322436 MeV |
| Maximum event-sum difference | 6.756817327868703e-13 MeV |
| PHITS histogram vs per-event HDF5 bin-count L1 | 0 / 150 bins |

HDF5 contains 1000 events, 123,322 particles, 5,977,725 tracks, and 5,978,912 deposits. It retains 4,797 exact-zero additions; there are no negative or invalid-energy rows. The scorer audit separately records one `NOT_SCORED_BY_PHITS` child candidate, zero unresolved, and zero unsupported additions. That candidate is not fabricated as a zero row.

HDF5 per-event total deposit mean is 1.059594627322437 MeV, sample SD 0.6249594592643488 MeV, standard error 0.019762953365425606 MeV. There are 87 zero-deposit histories and 275 histories at or above the tally upper edge 1.5 MeV. The comparison script reports overflow separately; those histories are not folded into the last tally bin.

Exact CLI and audit output are in the linked machine-readable artifacts listed above. No source Dumpall or input was modified.
