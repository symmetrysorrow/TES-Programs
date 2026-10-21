# PHITS 3.341陽電子scoreのソース追跡

監査したローカルPHITS 3.341ソースのファイルは `talls08.f`、`egs5edxde.f`、`egs5annih.f`、`egs2phits.f`。ソースそのものやバイナリはこの研究用ブランチへ含めない。

## 通常の陽電子輸送

`tdepstregEVENT` は輸送レコードの座標とcxyzから `tlngth` を作り、E/ECを取得する。ncol=9/11以外ではEGS5 `egs5edxde(ecc,tlngth,ein,mat,ityp)` を呼ぶ。ITYP=13の場合、`egs5edxde.f` は `ek=max(ein,emin(ityp))`、`lelke=eke1(m)*log(ek)+eke0(m)`、`dedx=pdedx1(lelke,m)*log(ek)+pdedx0(lelke,m)`、`eout=ein-dedx*step` を実行する。

この `[t-deposit] output=deposit` 条件では既存の`dexc_ene`が`OLDWT`とともに`trEVENT`へ加算される。Dumpallから求めたPDedx×endpoint chord×OLDWTは、同一seed、同一PHITS診断実行の1,795陽電子transport additionsと一致した。対象1,793件のraw-only独立計算では最大差1.39e-17 MeV。

## 消滅

EGS5の停止後消滅は2個の光子を作る。Dumpall側ではdead positron secondaryの`QCLUSTS(7)`は残余運動エネルギー、`QCLUSTS(5)`は静止質量であり、511 keV photon secondaryはITYP=14で別に輸送される。PHITSのtally dead charged particle分岐は子のQ7を`ein`、0を`ecc`として、条件を満たす子の局所scoreを生成する。通常枝ではQ7−0が加算値で、QCLUSTS(8)や静止質量は加算式に含まれない。

対象9個のdead positronのQ7と同じ実行のinternal tlvは全件完全一致した。これは消滅光子を同一点depositとして加える扱いを支持しない。二つの511 keV光子は輸送された光子であり、以後の光電吸収等で生じた荷電二次粒子を介したdepositだけが輸送過程上で計上される。

## 実機診断の健全性

10,000履歴の診断なしPHITS 3.341 runと、diagnostic `PRE/POST` traceを追加したisolated source buildを、同じseed・入力で比較した。両DumpallはSHA-256 `c11c51299c6c77a7e8a531a640ed72b31c8455fb2d0a17f2b03d3716a2815060`、両`[t-deposit]` histogramはSHA-256 `0d87d8b7041c54cf...`（元データとbaselineのファイル全体SHAはローカル比較ログ参照）で一致。輸送乱数系列・通常出力は診断追加で変化していない。

診断traceの内部数値は未公開ローカル検証出力に保管し、Gitへ追加しない。
