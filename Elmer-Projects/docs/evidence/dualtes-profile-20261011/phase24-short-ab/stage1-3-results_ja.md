# Phase24 SIGSEGV・Wall計測・数値比較の進捗

調査日: 2026-10-11
基準: `codex/dualtes-profile-evidence-20261011`、親コミット `f3e31e98`
対象: DualTES 4-step、MPI 6 ranks、MS-MPI / GNU Fortran Debug (`-O0 -g3 -fbacktrace -fno-omit-frame-pointer`)

## まとめ

1. Vector Assembly有効時のSIGSEGVは1 rankと6 ranksの両方で発生した。6 ranksのシンボル付きスタックでは、`CRS_PreparePositionCache()` と `CRS_GlueLocalMatrixCached()` のOpenMP `CRITICAL` から `libwinpthread-1.dll` に入った際に落ちた。キャッシュ位置検索やCSR値更新の行に到達する前の失敗である。これは発生箇所の特定であり、アプリ側のロック使用と実験ビルドのOpenMPランタイムのどちらが根本原因かはまだ切り分けていない。
2. Vector Assembly無効のDebug 6-rank測定では、`DefaultInitialize()` は各反復約0.003–0.005 sだった。一方、bulkと`DefaultDirichletBCs()`の担当時間はrank間で大きく異なる。全36 rank-iterationサンプルでbulkは4.460–15.019 s、Dirichletは0.008–9.862 s。rank 0/3はbulkが約4.5–5.5 sでDirichletが約9–10 s、rank 1/2/4/5はbulkが約12–15 sでDirichletは概ね0–1.7 sだった。全体assembly wallは13.727–15.040 sで、処理分担の不均衡が見える。
3. Vector ONは最初のassemblyで停止するため、温度、TES電流、RHS、行列値、収束履歴のOFF/ON比較は成立していない。速度向上も数値同等性も主張しない。75 ms本計算は再実行していない。

## 第1段階: 発生箇所

### 6 ranks

Debug symbols付き6-rank実行で、最初に捕捉したフレームは`CRSMatrix.F90:99`の`CRS_PreparePositionCache()`にある`!$OMP CRITICAL (PHASE24_POSITION_CACHE_SETUP)`だった。呼び出し元は`HeatSolve.F90`の`Phase24AssembleStaticElements()`で、並列要素ループの前にキャッシュ容量を準備する箇所である。

切り分けとしてこの準備ルーチンのcriticalだけをDebug用ソースで外すと、次は`CRS_GlueLocalMatrixCached()`の`CRSMatrix.F90:1174`にある同名criticalでSIGSEGVになった。呼び出し元は`Phase24AssembleStaticElements()`の要素ループ内だった。スタックはいずれもクリティカル区間入口のランタイム内で止まり、キャッシュbuild/hit後の数値更新には到達していない。

### 1 rank

未分割メッシュの1-rank SIFもVector Assemblyのstatic/prism fast assembly開始後にSIGSEGVで終了した。実行ログではrestart読み込みとメタデータ作成が完了している。ただし1-rank側では同じ箇所のシンボル付き有効スタックを採取できていないため、6-rankと同一の命令で落ちたとの断定は保留する。6分割メッシュをそのまま1 rankに渡した別試行は比較に使っていない。

### Debug bounds checkで見つかった別件

`-fcheck=all`を有効にした診断ビルドで、tetra metadata作成時に4要素の`nodeIds/eqIds`を6行の配列全体へ代入するshape mismatchを検出した。実験コピーでは`(1:4,i)`への代入に直すとこの検査エラーを通過した。この修正は今回のSIGSEGVの原因と確定していない。またruntime checksを有効にしたビルドはVector OFFでも通常経路と異なる早期失敗があったため、その実行結果は性能・再現性比較に用いていない。

## 第2段階: 同一Wall時計の区間測定

すべての記録区間は`RealTime()`を使い、rank-local CPU timeとは混ぜていない。MPI barrierは追加していない。rank値はstdoutへ直接出力し、各rank・各反復の値をログに残した。下表は36個（6 ranks × 6 nonlinear assembly calls）のサンプル統計で、平均はrank/反復を同じ重みで平均したもの。

| 区間 | 最小 (s) | 最大 (s) | 平均 (s) |
|---|---:|---:|---:|
| preinit (RadiationFactors/Tabulate等) | 0.002 | 0.004 | 0.003 |
| `DefaultInitialize()` | 0.003 | 0.005 | 0.004 |
| bulk前セットアップ | 0.000 | 0.001 | <0.001 |
| bulk assembly全体 | 4.460 | 15.019 | 10.599 |
| boundary assembly | 0.001 | 0.040 | 0.011 |
| finish処理 | 0.000 | 0.000 | 0.000 |
| `DefaultDirichletBCs()` | 0.008 | 9.862 | 3.631 |
| post/bookkeeping | 0.000 | 0.000 | 0.000 |
| 外側assembly全体 (同じRealTime時計) | 13.727 | 15.040 | 14.247 |

区間は直列の外側assembly内で測定した。preinit、initialize、prebulk、bulk、boundary、finish、Dirichlet、postの合計が外側assemblyに含まれる。ログ出力自体のコストを単独計測した対照試験は未実施であり、とくにDebug `-O0`の値を本番速度として扱わない。rankごとに同じ反復の外側assembly値がほぼ一致する一方、bulkとDirichletの負荷がrank間で入れ替わっているため、単純な要素数比較だけでは実行時間を説明できない。次はpartitionごとの境界条件数・自由度数と、Dirichlet内ループの要素数を採り、負荷不均衡の原因を確定する。

この計測は要求された大区間（初期化、bulk、boundary、finish、Dirichlet）を分けた。要素走査、材料物性評価、局所行列計算、CSR挿入、前処理生成の個別同時計測はまだ行っていない。

## 第3段階: OFF/ON比較と数値ゲート

| 条件 | 結果 |
|---|---|
| Vector OFF、DualTES 6 ranks | 4-step完走。Debug wall約91 s。数値出力あり |
| Vector ON、DualTES 6 ranks | 最初のassemblyでSIGSEGV。数値結果なし |
| Vector ON、DualTES 1 rank、未分割メッシュ | fast static/prism開始後にSIGSEGV。数値結果なし。該当箇所のシンボル付きスタック未採取 |

したがって温度場、TES平均温度・電流・抵抗・電力、RHS、行列sparsity/係数、線形残差、非線形履歴のOFF/ON差分は計算できていない。実行開始前後のstateファイルが同一hashのままなのは、初回組立てで停止し更新されなかったためで、数値一致の証拠ではない。fast pathは有効化しない。static matrix reuseもこの比較には含めていない。

## ビルドと適用範囲

Debug source/build/runtimeは`%TEMP%\phase24-debug-*`の隔離コピー。Debug時に追加したbounds修正、計測、critical切り分け用編集は一時ソースだけにあり、実行用の外部Elmer source treeおよび本番DLLを変更していない。実験用DLLのSHA256とビルド設定は`debug_manifest.json`に記録した。Vector OFF計測のSIFとrank別Wallログ、ON試行SIF、GDBスタックをこのディレクトリに保存する。

今回の実測で確定できるのは、Vector ONのSIGSEGV発生位置と、Vector OFF時の大区間rank別時間差までである。原因の根治、最適化後の速度と数値同等性は未解決。SIGSEGVを直し、同じ4-step入力でON/OFF両方が完走した後に、行列/RHS/場/回路量/収束履歴を比較してから改善効果を評価する。
