# DualTES Material 7 MATC評価の重複削減と数値検証

## 結論

Material 7の温度依存熱伝導率について、Generic assemblyでprismごとに行っていたMATC評価を、rank内のmaterial・unique nodeごとに評価する隔離実装を作成した。適用条件を満たさない場合は従来経路へ戻し、各非線形assemblyで現在の温度から値を再評価する。機能の既定値とVector AssemblyはOFF。

4-step、6-rank、同じDualTES restartでのDebug `-O0`単発比較では、rank 5 Bulk平均が13.473秒から3.759秒、外側wallが88秒から39秒に短縮した。OFF/ONの温度場最大絶対差は`1.3877787807814457e-16 K`、TES series/iteration CSVと6個のResult Normは一致した。初回非線形assemblyの行列・RHSキャプチャ38ファイルもSHA-256一致した。

これは隔離Debug環境での測定であり、Release性能、複数回ばらつき、4-step全assemblyの行列・RHS比較、HYPRE残差・反復数は未検証である。したがって本番統合・採用は未完了。パッチはレビュー用に保存し、本番DLL・本番設定は変更していない。

## 実装と適用条件

材料ごとにprism群のnode-to-slot対応表を作り、unique nodeの熱伝導率を材料単位で一括評価する。物性値は各非線形assemblyで現在の`Temperature`から更新する。データ構造はMPI rankごとに独立し、mesh変更時に解放・再構築する。

安全判定は限定的で、`Heat Conductivity`が`LIST_TYPE_VARIABLE_SCALAR_STR`、MATC依存名がTemperatureのみ、Lua/procedureなし、単一の熱伝導率定義の場合だけ対象にする。判定できない式はGeneric `ListGetRealArray`経路へfallbackする。行列、RHS、CSR挿入、時間積分、TES回路式は変更しない。

差分は外部Elmer診断ソースの一時コピーに適用した。リポジトリ内の製品ソースやDLLへは適用していない。差分の基点と適用範囲は`cache_patch_context_ja.md`を参照。

## 数値検証

同一MPI 6-rank分割、同じ4-step DualTES入力と同じrestart状態を用いた。Verify実行ではrank 1、2、4、5に属するMaterial 7 prism各2個について旧`ListGetRealArray`とキャッシュ評価値を直接比較し、差はすべて0。rank 0、3にはこのcaseのMaterial 7対象elementがない。

- 初回transient assemblyのA行列/RHSは全rankのbefore/afterファイルとサイズsidecar計38個のハッシュが一致。内訳は値ファイル24個とサイズメタデータ14個。
- TES 1/2のseries、iteration CSVはバイト一致。
- 全rankの4時刻温度場は最大絶対差`1.3878e-16 K`、最大相対差`8.227e-16`。前半の時刻は同一で、後半の差は浮動小数点最終bit水準。
- OFF/ONのHeatSolve Result Normは6値ずつ一致。
- HYPREログの線形反復数と最終残差は`-1`。実値を検証できていない。

大きなresult shardと行列/RHS payloadはコミットせず、rankごとのSHA-256と数値比較結果を`numeric_comparison.json`および`matrix_rhs_hashes.json`に保存した。生ログ、SIF、TES CSVは保持している。

## 性能結果 (Debug `-O0`)

OMP_NUM_THREADS=1、6 MPI ranks、OFFを先、ONを後に各1回実施した。

| 指標 | Cache OFF | Cache ON | 差 |
|---|---:|---:|---:|
| 外側wall (field output run) | 88 s | 39 s | 55.7%短縮 (2.26x) |
| rank 5 Bulk平均 (6 assembly計測) | 13.473 s | 3.759 s | 72.1%短縮 (3.58x) |
| rank 5 Material 7 material-eval区間 | 9.703 s | 0.0425 s | 99.56%短縮 |
| 全rank Assembly平均 | 13.491 s | 5.216 s | 61.3%短縮 |

rank 5はMaterial 7の36,023 prism/assemblyを持つ。旧経路では要素ごとに6 node値が要求され、理論上216,138 node-property値/assemblyに相当する。キャッシュでは22,024 unique nodeに整理し、Material 7の一括MATC評価呼出しは1回/assembly。構築時間は約0.038秒。別のVerify実行で計測したunique-node再評価区間は約0.98～1.15秒/assemblyで、追加コストは無視できないが、旧MATC区間より短い。

外側wallは単発測定で実行順も固定されており、OS cacheや起動状態の影響を分離できていない。Release性能への外挿はしない。

## ビルド失敗の切り分け

Ninja Fortran依存情報生成での失敗は、`.F90`の自由形式、`-cpp`前処理、CMake/Ninjaのmodule依存生成ではなく、Fortran frontend起動環境で再現した。`f951.exe --version`のWindows statusは`0xC0000139 (STATUS_ENTRYPOINT_NOT_FOUND)`。ucrt64/MSYS2のDLL directoryをPATH先頭へ置いた後、Vector Assembly OFF対照とキャッシュ試作の隔離Debug buildは成功した。保存ログとコマンドは`build_failure_diagnosis.log`、`build_off_control.log`、`build_trial_debug.log`にある。最初のNinja失敗時stderr全文は保存されていない。

## 未完了事項

1. Release最適化buildでOFF/ON性能を複数回交互に比較する。
2. 全4-stepの全非線形assemblyで行列・RHSを比較する。
3. HYPREの実線形反復数・最終残差を取得して収束を確認する。
4. single-pixel 4-stepで同じGeneric機構の数値・性能を検証する。
5. 温度以外に依存する材料式などの負例を確認する。
6. 上記ゲート後にソースを製品側へ統合するか判断する。

Vector AssemblyのSIGSEGVは独立した未解決問題。今回のパッチはVector Assembly OFFで動作する。75 ms解析は実施していない。

Release buildでの再検証と全timestep比較は[Release検証レポート](../matc-cache-release/report_ja.md)を参照。
