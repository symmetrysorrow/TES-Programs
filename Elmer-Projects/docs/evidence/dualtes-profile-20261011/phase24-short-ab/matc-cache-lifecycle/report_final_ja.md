# Phase24 MATC Node Cache 最終寿命・負例確認

## 判定

**結論は「実験候補」です。キャッシュ既定OFFを維持します。** 公式Elmer 26.1から01→02→03を再適用した隔離Releaseで、既存DualTES / single-pixelの4-step計算、AST単体試験、いくつかの実材料負例と数値一致を再確認しました。一方、同一プロセス内のModel/Mesh/material切替、solver再初期化、Cache OFF→ON→OFFは実行できていません。既存の時間依存式と6-rank構成での高速化は再現しましたが、未知構成に対する既定ONの安全ゲートは満たしません。

本作業では本番DLL、既存SIF、物理パラメータ、01/02/03パッチを変更せず、75 ms本計算およびVector Assemblyも実行・有効化していません。新しい試験は別のC:\\Temp buildと証拠用SIFで実施しました。

## ソース再現とビルド

- 元ソース: `ElmerCSC/elmerfem` `5a8de867068be0568f09af40fb90fee300dfbede`
- 適用順: `01_phase24_source_from_release_26_1.patch` → `02_matc_ast_purity_validator.patch` → `03_matc_cache_lifecycle.patch`
- 適用前blob、patch hash、適用後blob/hashの検査とRelease全ターゲットbuildは成功。再実行手順は既存の`replay_01_02_03_and_build.ps1`、今回のログは`final_replay_build.log`。
- 今回のRelease solver SHA256: `946ce238495957de1d95f6ca0e5f428ccaf966199df14cf7d83982ac22f40cae`
- HeatSolve.dll SHA256: `324c3168a16756f22ce1e80e85d0ac6939e2a8b751f651367f37636416ab8a12`
- libmatc.dll SHA256: `C745E793F5FFDEBD54404C6562C29C735F01BD470EE781EF75ACFCD6968BA2A7`
- 03はGeneric MATC node-slot map、node list、material ranges、数値bufferとカウンタから`SAVE`を外し、recursive `HeatSolver` activation内の状態としています。入口でclearし、mapは同activation内のassemblyでのみ使われます。物性値はassemblyごとに現在温度から再評価します。ローカルallocatableは正常なFortran procedure終了でスコープ解放されます。異常停止時のOSプロセス回収やライフサイクル変更時のruntime挙動は別の検証事項です。

## 寿命・メモリ検査

Debug隔離版を`-O0 -g -fcheck=all -finit-real=snan -finit-integer=-999 -fbacktrace -ffpe-trap=invalid,zero,overflow`で作り、6-rank / 4-step single-pixel restartをCache ON・Verify ONで実行しました。終了コードは0、MATC cached-vs-generic確認対象は全rankで最大絶対差0、配列範囲・浮動小数点runtime errorは記録されませんでした。4-step transientは同一solverプロセスで連続activationとassemblyを行い、各assemblyでの値再評価を確認しています。詳細は`debug_lifecycle_singlepixel.solver.log`とJSONにあります。

GCC AddressSanitizer capability probeはリンク時に`ld.exe: cannot find -lasan`となり利用できませんでした。Fortranのsentinel初期化は一般的な未初期化メモリ検出を保証しません。bounds-check実行も未解放メモリ、use-after-free、全経路での二重解放を網羅しません。ヒープリーク検出はSKIPです。

同一プロセスでのModel変更にはElmerの`ReloadInput`を使う案を確認しましたが、そのソース自身が再読込みでは全情報を読み込めず問題を導入し得ると注意しています。専用の安全な再初期化/Model変更ハーネスがない状態でこれを実行し、寿命PASSと見なすのは避けました。したがって下記は未実施です。

- solver再初期化、異なるModel、body-to-material割当て変更、MATC式変更
- node数不変のmesh connectivity変更、異materialがnodeを共有するfixture
- 同一プロセス内Cache OFF→ON→OFFとcache構築前early return

これらは`final_lifecycle_test_matrix.json`でSKIPまたはSTATIC_ONLYとしています。

## MATC安全性とfallback

Release版`libmatc`に対してAST単体テストを再実行し、Material 7の式、算術式、`abs(tx)`を受理し、`tx+ty`、座標/未登録symbol、未知関数、代入、制御構文、複数文、文字列、末尾トークン、不正括弧、空式、ユーザー関数を拒否しました。許可式は旧評価経路との数値一致を確認しました。全件`RESULT failures=0`です。

実SIFのMaterial 7を変更した隔離1-step入力では、次のfallbackをruntimeログで確認しました。

- `tx+ty`、`tx+x`、未登録symbol、未知関数、代入式、構文エラー: safety gateはcacheを拒否。Generic MATCも識別子/式エラーでMPI solverが非ゼロ終了しました。**cache fallbackは確認、物理計算成功はFAIL/該当せず**です。
- Material 7の入力依存を`Variable Time`と`Variable Coordinate 1`に変えた正しいSIF fixtureはcacheを拒否し、各1-step計算をGeneric経路で正常終了しました。
- Procedure conductivity、3成分異方性conductivity、複数conductivity成分: safe cacheを拒否し、Generic経路の1-step計算は正常終了しました。
- Lua conductivityの専用runtime fixtureは作れずSKIP。AST unitの禁止構文確認と同等とは扱いません。

MATC式のfail-closed動作および、Generic評価自体の不成立を別に記録しています。全実行のSIFと生ログは`negative_*`ファイルです。

## 数値回帰

同一Release DLL、同一6-rank、同じrestartと4-step SIFでDualTESとsingle-pixelをCache OFF/ON各1回実施しました。全実行のreturn codeは0、Result Norm、HYPRE反復回数/残差、TES series/iteration CSVのSHA256はOFF/ONで一致しました。

全6 rank・全node・全4時刻の温度比較:

| ケース | 最大絶対差 | bitwise差のあるnode値 | TES CSV |
|---|---:|---:|---|
| DualTES | 1.388e-16 K | 689 | byte-identical |
| single-pixel | 8.327e-17 K | 121 | byte-identical |

非常に小さい浮動小数点差は残りますが、保存結果として数値的に一致しています。比較の詳細、rank/timestep別node数、結果 shard hashは`final_release_temperature_comparison.json`です。今回の最終replayでは行列/RHS captureを再実施していません。前回03のmatrix/sparsity/RHS証拠はそのまま保持し、本回の全項目PASSには含めていません。

## Release性能

今回の再構成Release、Verify OFF、6 ranksで各ケース3回、OFF→ONを交互に実行しました。

| ケース | Cache OFF 平均±標本SD | Cache ON 平均±標本SD | 平均高速化 |
|---|---:|---:|---:|
| DualTES | 84.747 ± 0.903 s | 31.215 ± 0.217 s | 2.715x |
| single-pixel | 64.609 ± 0.119 s | 34.527 ± 0.201 s | 1.871x |

個別時間、result norms、HYPREデータ、DLL hashは`final_release_performance_runs.json`にあります。

旧03 Release DLLでも同じホスト、mesh/restart、rank、4-step input、Verify OFF、OFF→ON交互順で3回追加計測しました。旧03はOFF 86.337 ±1.480秒／ON 31.172 ±0.122秒、新版はOFF 84.747 ±0.903秒／ON 31.215 ±0.217秒です。ON平均差は約0.04秒、OFF平均差は約1.59秒で、いずれもこの小標本からキャッシュの寿命変更に起因する性能差とは断定できません。

以前の03実績89.481/32.383秒より、今回の旧DLL測定は両方速くなりました。旧DLLを同じ時期に測り直しても以前の値を再現しないため、旧測定と今回との差は少なくとも一部が時刻依存のmachine/runtime状態・外部負荷または計測条件差です。過去の同時刻CPU/OS負荷記録がないためその内訳は特定できません。02の74.985/26.333秒は別DLL hash/buildに基づくため、単純比較からmaterial cacheの有無やmap再構築を原因と推定しません。今回の同一DLL内の値は安定しています。

## 最終ゲート

- **実験候補: 到達。** 既存DualTES / single-pixelと、記載した負例セットでの動作・速度を確認しました。
- **安全なopt-in候補: 未到達。** solver再初期化と同一プロセスのModel/Mesh/material変更、OFF→ON→OFFを実行検証していません。
- **既定ON候補: 未到達。** Luaなど未試験の材料表現、モデル構成と再初期化経路も残っています。

従ってdefault OFFを維持します。次に必要な作業は、Elmerの実ライフサイクルを使うテスト専用 driver/harnessの用意、上記変更系を同一プロセスで実行した記録、Lua/未知材料fallback、debug heap/use-after-free検査、および必要なら最終行列/RHS回帰です。
