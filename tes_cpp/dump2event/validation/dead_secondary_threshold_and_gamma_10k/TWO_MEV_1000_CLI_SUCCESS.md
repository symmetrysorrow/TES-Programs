# 2 MeV / 1000履歴：通常CLI出力

状態: **完了、同一run内部加算と全件照合済み**。出力先はH:の新規専用ディレクトリ。入力Dumpallは読み取り専用。

CLIは通常の `dump2event.exe <dumpall> <h5> --physics pb-egs5 --hdf5-version 2 --physics-config <JSON> --stopping-table <CSV>`。実行開始は2026-10-21 21:46 JST。CLIログのelapsedは2349.8716 s。HDF5 v2を読み戻し、PHITS内部 `internal_additions.csv` と同一runのTALLYSEQを加算IDごとにjoinした。

PHITS内部の既知値は同一runの `internal_additions.csv` 5,978,912行。前段のraw-only照合は全加算数値を比較済みだったが、本検証ではscore gate適用後の通常CLI出力HDF5 v2を全件再照合した。

| 指標 | 実測値 |
|---|---:|
| PHITS内部加算 / HDF5 deposit | 5,978,912 / 5,978,912 |
| 同一加算IDの一致 | 5,978,912、欠落・重複0 |
| 内部にあってHDF5にない / 逆方向欠落 | 0 / 0 |
| イベント・粒子・領域ID不一致 | 0 / 0 / 0 |
| 有効座標比較 / 座標不一致 | 5,978,912 / 0、最大差0 cm |
| エネルギー最大絶対差 | 4.680023729664029e-14 MeV |
| エネルギー絶対差合計 / 1e-12超過 | 5.133758390738264e-10 MeV / 0行 |
| PRE増分対PHITS TLV最大差 | 1.1102230246251565e-16 MeV |
| HDF5総量 / PHITS内部総量 | 1059.594627322437 / 1059.594627322436 MeV |
| 履歴合計最大差 | 6.756817327868703e-13 MeV |
| zero scored addition / 負値 / invalid energy | 4,797 / 0 / 0 |
| スコアされない候補 / unresolved / unsupported | 1 / 0 / 0 |
| PHITS 150-bin histogram L1差 | 0 |
| 履歴deposit平均 ± 標準誤差 | 1.059594627322437 ± 0.019762953365425606 MeV |

`addition_id` は1から5,978,912まで連続。dataset件数はevents=1000, particles=123,322, tracks=5,977,725, deposits=5,978,912。履歴deposit分布はzero=87、overflow(1.5 MeV以上)=275。

完全join JSONは `H:\\TES-Programs-validation\\dead_secondary_threshold_and_gamma_10k\\two_mev_1000\\internal_hdf5_join.json`、HDF5件数監査とタリー比較は同じディレクトリの `hdf5_counts.json`, `hdf5_tally_comparison.json`, `hdf5_tally_comparison.csv`。Dumpall SHA256=`f84d0dd4624549bfcc1b54e705109aecc93c3b9a7ee345089d16d49baabc59dd`、HDF5 SHA256=`e90f120d548d3b128f6cac4b9068e8b749ffd013eb4997f51bd068bcf576d358`。
