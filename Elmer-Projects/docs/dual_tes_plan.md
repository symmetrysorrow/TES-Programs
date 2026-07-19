# 2TES・長尺吸収体構成(dual-TES)実装計画

日付: 2026-07-14 / 状態: **全フェーズ完了(A〜E)**。dual-TES 構成は生成・実行・物理検証済み。
20mm 吸収体スキャン構成(mesh_dual_20mm)追加済み

> **20mm スキャン構成 実施記録 (2026-07-15):**
> abs を 20mm・TES を両端(`tes_pitch = 19[mm]` → TES 中心 x=±9.5mm、abs 端から
> 0.25mm 内側)に。`build_mesh.py` にメッシュレシピ単位の `parameter_overrides`
> 機構を追加(1つの dual_tes ツリーをメッシュごとに異なる寸法で生成、ツリー複製不要)。
> グローバル `dabs_dx=20[mm]` / `tes_pitch=19[mm]` を既定にし、`mesh_dual_base` は
> `parameter_overrides` で 5mm/4mm に固定(旧構成を恒久再現)。新 `mesh_dual_20mm`
> (56,832ノード)。**発見・修正**: generate_project_geometry.py の内部 reconcile は
> 単数 `geometry` キーを再導出しないため、build_mesh 側で override 適用後に再 reconcile
> して解決済みツリーを注入(既存の mesh_overrides のみ機構にも潜在していた顕在化前のバグ)。
> 入射5点(左端から 30/60/90/120/150 目盛 ÷300、20mm基準 → 中心基準 x=-8/-6/-4/-2/0mm)
> を `case_dual20_pos30…pos150`(pulse、center明示x+y/z auto、restart from
> `case_dual20_steady`)として追加。`case_dual20_steady`: ITER=75 / NRM=0.15140160、
> L/R 対称性 1.21e-5。回帰: 単ピクセル SHA ビット一致、既存12ケースSIF無差分。
> 電流比較スクリプト `scripts/analysis/dual20_current_scan.py`(未実行ケースは警告スキップ)。

> **タイムステップ収束 実施記録 (2026-07-15〜16):**
> 5点スキャンの電流時系列が回復域で上下動(ギザギザ)。原因は**パルス後回復を
> 粗い1ms刻みで解いていた時間離散化誤差**(最大反復到達は0=カップリングは収束、
> よって反復不足でなく刻み起因)。回復域の刻みを段階的に細分化して収束を確認:
> pos30のL/R電流で successive-refinement RMS差(21–60ms)は
> 1ms→200µs=212/140nA、200µs→100µs=26/22nA、100µs→50µs=6/6nA、
> 50µs→10µs=2/1nA。**50µsで既に収束**(R離散ディップ深さ 0.687→0.525→0.441→
> 0.442→…µA、1msは過大)。ユーザ選択により最終は**入射直後のns段
> (1ns→10ns→100ns→1µs)を残し、以降40msまで10µs均一**(2112ステップ/ケース)で
> 全5位置を再実行(各~2.2h)。5ケース同一刻み=残差は共通モードで位置比較に無影響。
> 単ピクセル8ケースSIFはバイト不変。時系列比較スクリプト
> `scripts/analysis/dual20_current_timeseries.py`(L/R別2図、位置は青シーケンシャル)。
>
> **収束後の物理(粗い刻みの誤りを訂正)**: TES_L(近傍)は入射位置に強く依存
> (x=-8で深く速いディップ、中央で浅い)。TES_R(遠方)は**ディップ深さがほぼ
> 位置非依存**(全位置~190.35µA)で、位置情報は**onsetタイミング**に表れる
> (中央=最速着底、x=-8=最遅)。粗い1msで見えた「遠方ほど振幅大(反転)」は
> **主に時間離散化アーティファクト**で、収束解では消える。位置推定の頑健な指標は
> 振幅でなくタイミング、という結論はここでも成立。
>
> **ステージ構成の追加チューニング (2026-07-16〜17、pos30で反復検証)**:
> ユーザ要望で20-30ms域の段差(2µs→10µs→200µsの境界に起因)をさらに解消。
> 反復: 1µs/1ms → 2µs/200µs → **5µs(20.02-30ms均一)+100µs(30ms-)**。
> 中間段を撤廃し単一区間にしたことで20-30msの段差が完全に消失(ピーク精度は
> 5µsでも2µs比13nAのみの差、視覚上ほぼ同一)。**最終確定**: 5µs域は30msまで、
> 100µs域は**35msで打ち切り**(ユーザ指示、計算コスト削減。全5位置の応答は
> 20-35msに収まるため物理的情報は失わない)。5点スキャン全ケースをこの構成
> (2093ステップ、pos30のみ従来の700step版=100msまで温存)で再実行し、
> `dual20_current_L/R_timeseries.png` を最終確定。

> **リスタート回路状態の永続化 実施記録 (2026-07-15):**
> 過渡/パルスケースのリスタート由来過渡(定常場からリスタート時、UDF回路状態が
> `AverageTemperature=T0`(168.57mK)決め打ちで初期化され、実定常点(≈167.28mK、
> α大でこの1.3mKズレが初期電力を誤らせる)と不整合 → Aitken緩和が約20msかけて補正)を
> **定常回路状態の永続化**で解消。UDF に per-instance のオプション定数
> `KeyPrefix//'State File'`(`TES [L/R ]State File`)を追加: **定常モードは収束回路状態
> (AvgT/I/R/P/PrevI)を毎反復ファイルに上書き**(最後=収束値)、**過渡モードは初期化時に
> 存在すれば読んで種にし、書かない**(種を壊さない)。定数が無ければ Fatal せず T0
> フォールバック(=単ピクセル完全互換)。ビルダは **dual ケースのみ** に State File 定数を
> 出力(定常=書くパス、restartケース=restart_fromの定常が書く同一パス、メッシュdir内で
> .result と同じく run.py が移動しない場所)。**単ピクセル8ケースSIFはバイト不変**。
> 検証: 単ピクセル steady_3x = ITER 30 / NRM 0.15080609 厳密再現(後方互換)。
> case_dual20_pos30 のパルス前ベースラインスパンが **~64µK → 0.08µK(L)/0.10µK(R)**
> に平坦化、パルス応答は保存(L +42.7µK@20.22ms、R +18.0µK@26.02ms)。
> 実装ノート: 定常が書く PrevCurrent は T0初期化値のままだが、L/dt項(Dt=1msで
> L/Dt≈1.2e-5、バイアス項の約1000分の1)のため平坦性への影響は無視可(残留ドリフト
> 0.04µK)。単ピクセルへの同機構適用は将来課題(SIF再基準化が必要)。

> **フェーズD2 実施記録 (2026-07-15):**
> dual 3ケース(shunt_transient / pulse_center / pulse_offset)完走。物理検証:
> - **ゲート2(ベースライン平坦性)**: 整定後(t≥10ms)の温度スパン L=6.2µK /
>   R=2.9µK(<10µK)。初期の見かけ49µKはリスタート緩和2サンプル分で発振ではない。
> - **ゲート3(中央入射のL/R対称)**: 温度の最大相対差 5.1e-4(<1e-3)。
> - **ゲート4(偏心入射 x=+1.5mm)**: timing=PASS(R側が0.71ms先行、立ち上がり
>   0.09ms vs L側0.80ms)。amplitude は計画の素朴な予想と**逆**で、遠いL側の
>   ピーク(62.8µK)が近いR側(51.5µK)より約40%高い。
> - **VTU場データで独立確定**: CSV(UDFのAverageTemperature=ノード算術平均)は
>   VTU場と整定時サブµK・スパイク時≤3.5µKで一致。エネルギー保存は投入直後
>   E/E_pulse=1.000(4桁精度)。振幅逆転は実物理: 近接TESは低コンダクタンス
>   Stycast接触+高速局所シンクで浅いスパイクをすぐ排熱、遠方TESは遅延・拡散
>   平滑化された累積波面でより高く広いピーク(接触コンダクタンス低域通過+拡散平滑)。
>   **位置の頑健な指標はタイミング**であり振幅ではない。
> - `run.py` を L/R 系列CSV両対応に修正(NTFS大小無視で `_l_`/`_r_` として
>   書かれる問題も対処)。単ピクセルケースの挙動は不変。
> 検証スクリプト: `scripts/analysis/dual_series_analysis.py`(系列)、
> `scripts/analysis/dual_vtu_verify.py`(VTU場・エネルギー収支)。
> プロット: `generated/dual_baseline_series.png` / `dual_pulse_center_LR.png` /
> `dual_pulse_offset_LR.png`。
>
> **未了(将来課題)**: dual のベースメッシュ動作点(TES≈167.3mK)は 3x比 約1.1mK
> 低い(膜コンダクタンス離散化)。解像度収束を見るなら `mesh_dual_3x` が必要。
> L/R の R/I/P が定常でも1-2%ずれるのはメッシュノード非対称(TES_L 531 tet /
> TES_R 450 tet、体積は同一 4.00e-14 m³)由来のノイズ。

> **フェーズB 実施記録 (2026-07-14):**
> `scripts/support/mesh_names.py`(mesh.names パーサ)新設。build_cases.py の
> `BODIES` 定数と BC リテラル(1804/1104/1305/1204/1105/1205/1004)を撤去し、
> ケースのメッシュの mesh.names から名前引きで導出(材料ロールは `_L`/`_R` を
> 剥がした基底名、bath = `SiO2_2*__zmin` 全境界、モルタル3意味論ペア×サイド展開、
> 共有マスターは1BCに集約して複数スレーブが参照)。回帰: 既存8ケースSIF
> バイト一致。dual 導出: 19body・bath 2境界・モルタル6組を確認。
>
> **フェーズC 実施記録 (2026-07-14):**
> `tes_transient_heat_source.f90` を `TESCircuitModule`(`CircuitState` 型×3、
> 全16状態フィールド)+共通 `TESCircuitCompute`(数式・実行順序不変)に
> リファクタ。エントリポイントはモジュール外 plain FUNCTION 3本:
> `TESTransientHeatSource`(無プレフィックス、unit 91)/ `...L`(`TES L `、92)/
> `...R`(`TES R `、93)。回帰: steady_3x = ITER 30 / NRM 0.15080609 を厳密再現。
>
> **フェーズD1 実施記録 (2026-07-14〜15):**
> ビルダ dual 対応(TES_L/R → BF1/2 に `TESTransientHeatSourceL/R`、Constants は
> `TES L/R ...` フルセット、Series File は `_L`/`_R` 挿入)。パルス中心スキーマ
> 拡張(`"center": {"x": "1.5[mm]", "y": "auto", "z": "auto"}`、auto成分=abs重心、
> Discrete Norm は解決後中心で自動計算)。dual 4ケース追加・生成、既存8ケースは
> バイト一致維持。`case_dual_steady`: ITER=67 / NRM=0.15042566 で収束、
> TES_L=0.1672581 K / TES_R=0.1673430 K(相対差 5.07e-4 <1e-3 ゲート通過)、
> 単ピクセル同解像度の一時実験値 0.1673230 K を挟む正しい転移動作点。
> ベースメッシュの動作点は 3x比 約1.1mK 低い(解像度効果)。
> **教訓**: ASCII `.result` の値は**ノード昇順**(perm逆引きではない)。
> 誤読すると正しい解が「冷分岐転落+浮遊ノード」に見える(誤診断1往復の実費)。
> メッシュ健全性は連結成分解析で確認済み(dual 5成分=単ピクセル3成分の相似形、
> 浮遊クラスタ・重複ノードなし)。

> **フェーズA 実施記録 (2026-07-14):**
> `geometries` レジストリ化(single_pixel 無変更移動 + dual_tes 新規19body)、
> `meshes.*.geometry` 参照、`mesh_dual_base`(54,156ノード)生成。
> `generate_project_geometry.py` の単スタック前提(body名・タグ101–109、
> membrane分割、モルタルretag)を sides ループに一般化
> (`tag = 101 + 9*side + role`、単ピクセルは従来値を厳密再現)。
> 回帰: 単ピクセル `gmsh/project.msh` SHA256 ビット一致
> (b42c9c3f…)、既存8ケースSIF無差分。
>
> **重要な発見**: 旧 `fragment_mortar_interfaces` は**最初から no-op だった**
> (ツールシリンダーが literal y=0 生成でジオメトリは y=1mm、かつメッシュ生成後に
> 実行)。既存メッシュの `abs__zmin_free`(全面重複)/`TES__zmax_free`(側壁
> スリバー)は偽物。物理的な接触制限は「スレーブ面=Stycast 真円ディスク」で
> 担保されていた。dual では意図された機能を pre-mesh で正しく実装:
> 0.16µm TES 膜は OCC ブーリアン不能のため footprint を 2D 分割してから一括
> extrude(コア+リング共形、体積和厳密)、abs 底面は 2D ディスクインプリント
> (シーム 0.082 rad 回転で ElmerGrid merge の円周ノード直結を回避)。
> 接触パッチ8面すべて厳密に π·(249µm)² = 1.9478e-7 m²。`abs__zmin` は
> 両側の接触円2枚を単一グループとし、L/R スレーブが同一マスターBCを参照する
> 設計。単ピクセルは SHA ゲート凍結のため旧経路を温存(偽 `_free` は
> データ駆動化で参照しないこと)。

## 0. 目的と決定事項

単ピクセル構成(abs が TES 基板上に載る)に加えて、**1本の細長い abs の両端に
TES 基板スタックを置く 2TES 構成**を追加する。電気回路(バイアス/シャント/
インダクタンス)は左右で独立に 2 系統持つ。

ユーザ決定事項:

- **既存単ピクセル構成と共存させる**(ジオメトリをレジストリ化。既存 8 ケースは
  回帰テストとして機能し続ける)。
- **寸法は仮置きで開始**: abs = 5×1×0.7 mm、TES ピッチ(左右 TES 中心間)= 4 mm。
  全て `parameter_expressions` でパラメータ化し、後から JSON 編集+メッシュ再生成
  のみで変更可能にする。基板チップ(3×6 mm)・膜・TES・Stycast は現行と同一断面。

## 1. 全体像

変更は 4 層に分かれる。フェーズ A〜D を順に実施し、各フェーズ末尾の回帰ゲートを
通過してから次へ進む。

| フェーズ | 内容 | 回帰ゲート |
|---|---|---|
| A | ジオメトリレジストリ化 + dual_tes ジオメトリ + メッシュ生成 | 単ピクセル `gmsh/project.msh` SHA256 ビット一致 + 生成 SIF 無差分 |
| B | build_cases.py の body/BC データ駆動化(mesh.names 由来) | 既存 8 SIF バイト一致 |
| C | UDF の回路多インスタンス化 + DLL 再ビルド | steady_3x = 30 反復 / NRM 0.15080609 再現 |
| D | dual ケース定義・生成・実行・物理検証 | 対称性(中央パルスで L/R 一致)ほか下記 |
| E | ドキュメント整備 | — |

実装は Sonnet サブエージェントに委譲し、フェーズ間のレビューと統合判断は
指揮側(メイン会話)が行う。

## 2. フェーズ A: ジオメトリレジストリと dual_tes ジオメトリ

### A-1. スキーマ変更(elmer_project.json)

- `geometry`(単一)→ `geometries` レジストリへ:
  `geometries.single_pixel` = 現行の geometry ツリー(**内容は無変更**)、
  `geometries.dual_tes` = 新規。
- `meshes.*` に `"geometry": "<name>"` を追加(既存 2 メッシュは
  `single_pixel`)。新メッシュ `mesh_dual_base` を追加(`dual_tes`、
  ベース解像度、ElmerGrid 引数は既存と同型で `-out mesh_dual_base`)。
- 新パラメータ(命名は既存流儀に合わせる):
  `dabs_dx = 5[mm]`(長尺 abs)、`tes_pitch = 4[mm]`。y/z 方向・各層厚は
  既存パラメータを再利用。

### A-2. dual_tes ジオメトリ内容

- `abs`: 1 個、中心 x=0、`dx_expr = dabs_dx`(y=1mm、z は現行式と同じ)。
- 左右スタック(サフィックス `_L` / `_R`、`x_expr = -tes_pitch/2` / `+tes_pitch/2`):
  TES、Stycast、Membrane、SiO2_1、Si_1(+sub)、SiNx(+sub)、Si_2(+sub)、
  SiO2_2(+sub)。ブーリアン構造(membrane くり抜き)は現行の複製。
- 期待される body 数: 1 + 2×9 = 19。境界 ID は既存規約
  (body=100+index、boundary=1000+100×index+面)が自動で振られる。
- 接触(モルタル)面: vendored contact_detection が自動検出。左右それぞれに
  TES↔Membrane_SiNx、Stycast↔TES、abs↔Stycast の 3 ペアが立つこと。

### A-3. コード変更

- `build_mesh.py`: レシピの `geometry` 名でレジストリから選択し、解決済み JSON
  (`generated/_mesh_build_input.json`)に従来キー `geometry` として注入
  → `generate_project_geometry.py` と vendored loader は無変更で済む想定。
- `reconcile_project.py` / `sync_elmer_parameters.py` / ドリフト検査:
  `geometries` 配下の全ツリーの `*_expr` を評価するように更新。

### A-4. 回帰ゲート(A)

1. 単ピクセルメッシュ再生成で `gmsh/project.msh` の SHA256 が
   `b42c9c3fbbeaeeb46d2971b15268bb694f39492d040187ddb049873cfd9b7bea` に一致。
2. `python sync_elmer_parameters.py` 後、`generated/cases/` に git 差分なし。
3. `python build_mesh.py mesh_dual_base` が成功し、`mesh.names` に 19 body と
   左右分のモルタル面が揃う。断面プレビュー(xz/yz)を出力して形状を目視確認。

## 3. フェーズ B: body/BC のデータ駆動化

現在 `build_cases.py` は `BODIES`(100〜109→材料)と BC の target ID
(1804/1104/1305/1204/1105/1205/1004)をハードコードしている。これをケースの
メッシュの `mesh.names` から導出する形に変える。

- `mesh.names` をパースして body 名→ID、境界名→ID を得る。
- 材料ロールは body 名からサフィックス `_L`/`_R` を剥がした基底名で引く
  (abs→Pb、TES→TES、Stycast→Stycast、SiO2_*→SiO2、Si_*→Si、SiNx→SiNx、
  Membrane_*→Membrane)。
- モルタル BC は意味論ペアの固定リスト
  (TES__zmin↔Membrane_SiNx__zmax、Stycast__zmin↔TES__zmax、
  abs__zmin↔Stycast__zmax)を、存在するサフィックスごとに展開。
  bath BC は `SiO2_2*__zmin` の全境界を対象にする。
- **出力順・整形は現行 SIF と完全一致させる**(単ピクセルではバイト一致が
  ゲート)。

### 回帰ゲート(B)

- `python sync_elmer_parameters.py` 後、既存 8 ケースの SIF が git 差分なし
  (バイト一致)。

## 4. フェーズ C: UDF の回路多インスタンス化

`tes_transient_heat_source.f90` の回路状態(SAVE スカラー群)を派生型
`CircuitState` の配列(3 要素)に括り出し、共通サブルーチンに
インスタンス番号を渡す。

- 公開関数:
  - `TESTransientHeatSource` — インスタンス 1、定数プレフィックスなし
    (**現行と完全互換**。単ピクセルケースは無変更で動く)。
  - `TESTransientHeatSourceL` / `TESTransientHeatSourceR` — インスタンス 2/3、
    定数は `TES L ...` / `TES R ...`(例: `TES L Bias Current`)、系列 CSV は
    `TES L Series File` / `TES R Series File`。
- 状態には LastTimeStep/LastNonlinIter/sweep 集計/Aitken(ω, cap, 残差)/
  FileStarted まで**全て**含める(片側のコミットが他方を踏まない)。
  ロジック本体(暗黙結合+Aitken 緩和)は一切変更しない。
- `AbsorberWindowPulseHeatSource` は無状態なので変更不要。
- 再ビルド: `elmerf90 tes_transient_heat_source.f90 -o tes_transient_heat_source_t0.dll`
  (Elmer 26.1、PATH は `C:\Program Files\Elmer 26.1-Release`)。

### 回帰ゲート(C)

- `python run.py case_tes_steady_3x_refined` を再実行し、収束反復数 30・
  NRM 0.15080609 を再現(リファクタが数値経路を変えていない証明)。

## 5. フェーズ D: dual ケースの定義・生成・実行

### D-1. ビルダ拡張

- メッシュに `TES_L`/`TES_R` body があるケースでは Body Force を 2 本
  (L→`TESTransientHeatSourceL`、R→`...R`)+パルス用 1 本の構成で生成。
  Constants には L/R 両プレフィックスのフルセットを出力(初期値は共通の
  回路パラメータ。将来ケース spec で片側上書き可能な構造にしておく)。
- パルス中心のスキーマ拡張: `"center": "auto"`(abs 重心)に加えて
  `"center": {"x": "1.5[mm]", "y": "auto", "z": "auto"}` 形式の明示指定
  (auto 成分は重心座標)。`Pulse Discrete Norm` の自動計算は中心指定に追従。

### D-2. 新ケース(mesh_dual_base、タイムステップ段構成は既存 pulse 系を流用)

| ケース | 内容 |
|---|---|
| `case_dual_steady` | 定常。restart 起点 |
| `case_dual_shunt_transient` | 無パルス過渡(ベースライン平坦性) |
| `case_dual_pulse_center` | 中央入射(x=0)。L/R 応答一致の対称性検証 |
| `case_dual_pulse_offset` | 偏心入射(例 x=+1.5mm)。L/R 非対称の位置感応検証 |

### D-3. 検証ゲート(D)

1. `case_dual_steady` 収束、左右 TES の平均温度が一致(対称性)。
2. `case_dual_shunt_transient` のベースラインが平坦(既存基準: span 数 µK)。
3. `case_dual_pulse_center` で L/R 系列 CSV が数値誤差内で一致。
4. `case_dual_pulse_offset` で R 側の応答が先行・増大し L 側が遅延・減少
   (定性)。エネルギー保存(入射 1332 keV に対する両 TES+熱浴への収支)は
   既存の解析スクリプトがあれば流用。

## 6. フェーズ E: ドキュメント

- `README.md` / `docs/redesign_plan.md` に dual-TES 構成の節を追記。
- 本計画書に実施記録を追記(redesign_plan.md と同じ流儀)。

## 7. リスク・注意点

- **UDF リファクタの最大の罠**: タイムステップ/非線形反復の遷移検出を
  インスタンスごとに独立させること。共有すると片側の状態コミットが欠落する。
- 定常の収束反復数・NRM は dual 構成では当然変わる。新基準はフェーズ D で採取。
- dual メッシュはベース解像度でも体積が約 2 倍。パルスケースの実行時間は
  数十分オーダーを見込む(バックグラウンド実行)。
- 左右チップ間ギャップ: ピッチ 4 mm・チップ幅 3 mm で内側エッジ間 1 mm。
  abs(5 mm)は両チップに 0.5 mm ずつ張り出す。Stycast 柱はチップ中心
  (=TES 中心)に置く。寸法変更時は干渉(チップ同士の重なり)に注意。
- `mesh.names` の面番号規約(…04=zmin、…05=zmax)は vendored ビルダ由来。
  データ駆動化はこの規約でなく**名前**(`TES_L__zmin` 等)で引くこと。
