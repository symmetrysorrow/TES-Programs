# 10,000履歴の陽電子未解決とヒストグラム差

## 判定

この独立seedの150-bin histogram差は、PHITS 3.341の内部 `trEVENT` 加算を履歴IDと加算順で記録した診断runとの照合により、HDF5 v2に残っていた9履歴の陽電子寄与欠落で説明できた。13-bin、L1=18は、その診断runから得た履歴別depositを同じ半開区間 `[eb(i), eb(i+1))` で再binしたときに再現する。HDF5の342 UNRESOLVED行だけが欠落ではなく、陽電子EGS5輸送1,460行が `E==EC` を理由に誤って既知ゼロとしていた。

## 排他的な行分類

| 対象 | 件数 | 生Dumpallに必要な情報 | 判定 |
|---|---:|---|---|
| 既存HDF5 `unsupported_primary`：陽電子NCOL=14/JCOLL=14 | 326 | E, OLDWT, endpoint chord, PB medium 2 PDedx table | 全件再現 |
| 既存HDF5 `unsupported_primary`：陽電子NCOL=13/JCOLL=2 | 7 | 同上 | 全件再現 |
| 既存HDF5 `unsupported_secondary`：dead陽電子 | 9 | JCLUSTS(4), QCLUSTS(7), child position | 全件再現 |
| 既知ゼロだった陽電子輸送（実際の `tlv>0`） | 1,460 | 同上のEGS5陽電子輸送情報 | 全件再現 |
| 新たに対象にした陽電子寄与 | **1,802** | — | 排他的区分の合計 |

最初の3区分が従来の342 UNRESOLVED行（333 primary + 9 secondary）である。1,460件は旧HDF5で energy-valid なゼロ行だったため342件には含まれず、別に数えた。PHITS internal traceと順序joinした全1,954,742加算のうち、対象陽電子の4区分は計1,802件。

未解決行があった9履歴はNOCAS `507, 915, 922, 2085, 2291, 6121, 6936, 8209, 9571`。履歴別の事前既知小計、internal合計、修復後HDF5合計は [`UNRESOLVED_HISTORY_IDS.csv`](UNRESOLVED_HISTORY_IDS.csv) に示す。

## エネルギー内訳

- 1,460件の `E==EC` 陽電子輸送: raw-only予測合計 0.6478282322839034 MeV。
- 333件の陽電子primary輸送: raw-only予測合計 0.07078472447285694 MeV。
- 9件のdead陽電子secondary: QCLUSTS(7)合計 0.008936438570298348 MeV。
- 新たに回収する合計: **0.7275493953270587 MeV**。これは旧HDF5の9不完全履歴で内部traceとの差として観測した合計 0.727549395327054 MeVとbinary64集積精度内で一致する。

PDedx再構成は1,793 primary record全件で同一計算のPHITS internal `tlv` と比較した。許容差は事前固定の1e-12 MeV。全件pass、最大絶対誤差1.3877787807814457e-17 MeV。dead陽電子9件はQCLUSTS(7)とinternal `tlv` の差0件、最大誤差0 MeV。

## 物理的根拠

PHITS 3.341 `talls08.f:tdepstregEVENT` は粒子ITP=12/13/14に `egs5edxde(ecc, tlngth, ein, mat, ityp)` を呼ぶ。`egs5edxde.f` は陽電子ITP=13なら `pdedx1(index)*ln(E)+pdedx0(index)` を計算し、`eout=ein-dedx*step` とする。Deposit event score側の `dexc_ene` は旧重み `OLDWT` を乗じて `trEVENT` に加算される。Dumpallのraw `E`、endpoint chord、OLDWT、material/medium mapping、およびHATCHの明示的PDedx係数から同式を再計算できる。

Dead charged-secondary分岐は `JCLUSTS(4)<0`, `JCLUSTS(5)!=0`, `QCLUSTS(7)>eletb(1)` を判定し、`ein=QCLUSTS(7)`, `ecc=0`, 通常枝では `tlv=ein-ecc` とする。QCLUSTS(8)とOLDWTはこの枝の加算式で掛けられない（ソースコメントで重み乗算は削除済み）。対象9件はこのスコアしきい値より十分大きい。QCLUSTS(7)は運動エネルギーであり、陽電子静止質量や消滅光子511 keVをこの点depositへ加えていない。

## 範囲

検証済み範囲はPHITS 3.341、Pb region 1000、Dumpall material ID 2 / EGS5 medium 2、1 keVの電子・陽電子cutoff、単位粒子重みの同一10,000履歴seedである。重みが履歴中に変わる場合、別材料、他region、別のscore unit、ユーザー定義deposit関数への一般化は未検証。raw/source traceと係数表の全行は `private/` に置き、公開用報告には含めない。

## Exact bins changed by the prior partial HDF5

The prior partial histogram had 13 signed bin differences (bin-index: reconstructed minus PHITS): `11:-1, 12:-1, 14:-2, 15:+3, 18:-1, 27:-1, 30:+3, 45:-1, 55:+1, 69:-1, 70:-1, 81:+1, 82:+1`; total L1 count difference 18. The new internal PRE trace allowed all 9 event histories and every tally addition to be identified directly; no matching was inferred from histogram shape. After reconstructing the 1,460 false-zero positron tracks and 342 unknown rows, all 150 HDF5-vs-PHITS bin counts match.
