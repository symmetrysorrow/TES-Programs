# 全加算・履歴・位置の検証結果

## 実行

- PHITS: 3.341 sourceを隔離した診断buildで同一seedを再実行。実行日時 2026-10-21 11:57:40–12:05:25 JST、exit code 0、`END` markerあり。
- 元10,000履歴Dumpall: `dump2event_physics_validation/overnight_physics_reconstruction/runs/n10000_independent_seed/dumpall.dat`。2,783,824,000 bytes、SHA-256 `d254b3afd13230be84ad8f733c8aaa89811724f832a506a0d5b421d051c06b1c`。修正前後でサイズ・更新日時に変化なし。
- PHITS診断の再生ではなく、10,000履歴の元通常計算と同一seed・同一設定の再計算をした。計算後Dumpallとbaseline Dumpallは同じSHA-256 `c11c51299c6c77a7e8a531a640ed72b31c8455fb2d0a17f2b03d3716a2815060`。`[t-deposit]`出力もbaselineとdiagnosticで同じSHA-256 `0d87d8b7041c54cfd10e610ac0653ff1eda1844a31ca0691126ce82ff0c5f110`。
- internal tally trace: 1,954,742 score additions、10,000 histories。traceはこの1 runの`PRE`値と区間端点を順に記録する診断出力。

## 同一実行の加算照合

生Dumpallから独立に計算した実験用ledgerと同一計算のPHITS internal additionsを、順序、NOCAS、particle、regionで照合した。

| 指標 | 結果 |
|---|---:|
| Dumpall structural records / histories / terminal record | 1,983,461 / 10,000 / NCOL=2 |
| 加算ID順のPHITS internal score calls | 1,954,742 |
| raw-only reconstructed additions | 1,954,742 |
| UNRESOLVED / negative | 0 / 0 |
| energy mismatch > 1e-12 MeV | 0 |
| region mismatch | 0 |
| positive-deposit location mismatch > 1e-12 cm | 0 |
| maximum per-addition absolute energy difference | 5.865533753146579e-15 MeV |
| maximum checked position difference | 0 cm |
| positive-energy additions | 1,943,292 |

PHITS internal per-history sumとHDF5 v2の履歴別合計は、全10,000履歴でbin assignment一致、最大合計差5.340172748447003e-14 MeV。合計はinternal 332.78211273165505 MeV、HDF5 332.7821127316546 MeV。

## Histogram comparison

同じPHITS run、same 10,000 history set、region 1000、unit=3 (MeV/source)、0–1.5 MeV、150 bins、0.01 MeV幅を比較した。positive in-range historiesは双方532。HDF5履歴和をPHITSの同じ半開区間へ入れたとき、bin count差は0、L1=0。PHITSテキストの表示値（4桁丸め）との浮動小数点差はL1 `1.554312234475219e-14`。いずれのビンもcount差なし。最大binは1.17–1.18 MeVに36履歴、nonzero binsは121。

修復後HDF5 v2は587,741,058 bytes、SHA-256 `82fd2704c670d383cae81b66ff50c5aa30b27e3d599e8ef50bcebbadd49535ea`。

旧HDF5 v2は342 unresolved行を持ち、positive event数は偶然532でも、13 bins/L1=18が異なった。補正HDF5 v2は全加算を同一内部traceで照合し、PHITS histogramにも一致した。比較は過去の別seed histogramと混ぜていない。

## 位置の意味

1,943,292 positive score rowsのDumpall輸送区間のstart/endはPHITS `PRE` traceのraw区間端点と全件一致した。dead secondary 1,319点はDumpall QCLUSTS位置を使い、該当内部score regionはすべて1000、診断区間終点と一致した。PHITS診断hookは加算ごとの変換後代表点を別出力していないため、dead-secondary pointの独立座標検証はraw child coordinate + `trnsxx/tregck`のソース分岐 + internal region/endpointによる確認である。連続区間のPHITS代表点を作り出してはいない。

**判定:** この10,000履歴のunit-weight Pb region 1000 sliceは加算値・履歴・regionが一致。その他のseed/重み/材料/generalized voxel transportへの合格を意味しない。
