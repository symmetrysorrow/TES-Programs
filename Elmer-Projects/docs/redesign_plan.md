# プロジェクト構成の再設計案

日付: 2026-07-14 / 状態: フェーズ0・1・2・3 実装済み(下記)

> **フェーズ3 実施記録 (2026-07-14) — 外部GUI(Thermal-and-Electoric-Sim)を完全に無視して良い
> と確認されたため実施:**
> gmshジオメトリスタックをベンダリング(`scripts/support/vendored/geometry/`:
> gmsh_builder / contact_detection / primitives / spec / body_semantics は verbatim、
> templates は freeform_boxes サブセットのみ)。4箇所のモンキーパッチは
> `TesGmshBuilder(GmshApiBuilder)` の通常のオーバーライドに統合。GUI形式の
> 全機能パーサ io.py(1,108行)は本プロジェクト用の軽量ローダ
> `vendored/geometry/loader.py`(~110行)で置換。**これで外部リポジトリ依存はゼロ**。
> 等価性検証: 外部経路は決定的(再生成で bit 一致)であることを確認した上で、
> ベンダリング経路・式のみ化スキーマの双方で `gmsh/project.msh` の SHA256 が基準
> (B42C9C3F…)と**ビット単位一致**。
> JSONスキーマ最終整理: `parameters`(数値)/`parameter_units`/`role_map`/
> `material_assignments` と他シミュレータ用セクション(boundaries/probes/sources/
> circuit_layout/resistance/solver/contacts/analysis_planes)を削除、materials は
> expression のみ、geometry は `*_expr` のみ(数値は reconcile がビルド時に導出)。
> 28.6KB→18.4KB。ドリフト検査は「双方に存在して不一致」のみ報告に変更
> (欠落=設計どおり導出専用)。ケースSIF生成物はバイト単位で不変を確認。

> **フェーズ2 実施記録 (2026-07-14):**
> `run.py`(restartチェーン自動解決、実行、`results/<case>/` への出力整理、
> 入力ハッシュ付き `manifest.json`)と `build_mesh.py`(`meshes` レジストリ →
> ジオメトリ生成+ElmerGrid+`PROVENANCE.json`)を追加。`meshes` レジストリを
> `elmer_project.json` に追加し、ケースビルダがメッシュ存在を検証。既存2メッシュには
> `--record-only` でPROVENANCEを記録(3xのレシピは事後再構成のため `recipe_verified:
> false`)。**計画からの変更**: gmshビルダ群のベンダリングは見送り(外部リポジトリで
> ~3,000行+推移的依存があり、GUI側で開発継続中のコードの複製は保守コスト過大)。
> ジオメトリ生成のみ明示的な外部依存として残す。`Results Directory` キーワードは
> 使わず、run.py が実行後に出力を移動する方式(restart interfaceの `.result` は
> メッシュディレクトリに残す)。`role_map`/`elmer_overrides` はジオメトリ生成が
> 実使用しているため分離しない。

> **フェーズ1 実施記録 (2026-07-14):**
> 全8ケースを `elmer_project.json` の `cases` セクション(template: steady/transient/pulse)
> に移行し、`scripts/support/build_cases.py` が自己完結SIFを `generated/cases/` に生成。
> 手書きルートSIFと旧 `tes_case_*.sif` フラグメントは削除。パルス中心(abs重心)と
> `Pulse Discrete Norm` は `scripts/support/mesh_quantities.py` がビルド時に自動計算。
> UDFは定数欠落で `Fatal`(暗黙デフォルト全廃)、旧 `AbsorberPulseHeatSource` 削除。
> ケース値は式(`"20.02[ms]"` 等)で記述可、エネルギーはkeV基底→J変換。
> **回帰で判明した旧チェーンの欠陥2件**: (1) `Real $ I0` のMATC文字列置換が
> I_0 を~7桁に切り詰めていた(生成SIFは完全精度; 平衡ノルム差 1.4e-7 相対)。
> (2) 旧 `AbsorberPulseHeatSource` の `Y0=0` ハードコードは現メッシュのabs中心
> (y=1mm)を外していた(自動重心で恒久解消)。回帰基準: steady_3x = 30反復 /
> NRM 0.15080609。SIF数値の型付け注意: 裸の整数はinteger型エントリになるため、
> ソルバ実数スカラーは小数点付きで出力(`fmt_real`)。

> **フェーズ0 実施記録 (2026-07-14):**
> 式評価器を `scripts/support/vendored/` にベンダリングし、`reconcile_project.py` と
> `sync_elmer_parameters.py` から外部リポジトリ依存を除去。JSONの自己書き換えと
> `autosave/` を廃止(ずれは警告+`--fix-numerics` での明示更新に変更)。
> `generate_project_geometry.py` は解決済みJSONを `generated/_project_resolved.json`
> に書いてローダーへ渡す方式に変更(主ファイル不変)。ただし gmsh ビルダの
> 外部依存はフェーズ2まで残る。`TES_volume` は `TES_Au_dx*TES_Au_dy*TES_dz` から
> 導出、`L_tes = 12.3[nH]` をパラメータ化し `TES Inductance` として全transient系
> SIFのConstantsへ(UDFのハードコード撤去)。回帰: 定常3xケースで確認。

## 1. 現状の問題点(診断)

### P1. 外部リポジトリへの絶対パス依存

`scripts/support/reconcile_project.py` と `generate_project_geometry.py` は
`D:\Github\Thermal-and-Electoric-Sim` を `sys.path` にハードコードで注入し、
`core.project.dimensioned_expression`(単位付き式の評価)や
`core.geometry.gmsh_builder` をインポートしている。

- `elmer_project.json` の意味論(`"715[uA]"` が何を意味するか)が**外部リポジトリの
  Pythonコードで定義**されており、このリポジトリ単体では値を確定できない。
- 外部リポジトリの改変・移動・削除で、ここのビルドが黙って壊れる。

### P2. 同じ値の三重持ちと自己書き換え

`elmer_project.json` は各値を3系統で持つ:

- `parameters`(SI数値)+ `parameter_units` + `parameter_expressions`(式)
- `materials.*.{nominal, expression}`
- `geometry` の `x` / `x_expr` ペア

`sync_elmer_parameters.py` はロード時に式から数値を再導出し、**差があればJSONを
書き戻す**(自己書き換え + `autosave/` スナップショット)。真実の所在が
「式」なのか「数値」なのか構造上は判別できず、Pythonを通さない編集は不整合を生む。
git差分もノイズ化する。

### P3. ケース定義がPythonにハードコード、主力ケースがJSON管理外

- `sync_elmer_parameters.py` の `elif` 連鎖は4ケース
  (`shunt_internal / shunt_transient / pulse_1332kev / constant_power`)のみ。
- 現在の主力である **3x 系ケース(steady / transient / pulse_20ms)は完全手書きSIF** で、
  以下がJSONの管理外に散在している:
  - パルス定数(エネルギー・時刻・幅・σ・中心座標)
  - タイムステップ多段構成(10ステージ配列)
  - solver設定のケース別上書き(例: 非線形反復 25 は共有値 15 の手書き上書き)
  - `Pulse Discrete Norm`(**メッシュ依存量をオフライン手計算**して貼り付け)
  - `TES Series File` などの出力先

### P4. UDFのデフォルト値による二重管理

`tes_transient_heat_source.f90` は全定数に `GetConstReal(..., Found)` +
ハードコードのフォールバックを持つ。

- SIFで定数が欠落しても**古い既定値で黙って動く**(バグの温床)。
- インダクタンス `L = 1.23e-8` は **UDF内ハードコードのみ**(SIF/Constantsに存在しない)。
- `TES_volume = 4.0e-14` は `sync_elmer_parameters.py` 内のリテラル
  (本来 `TES_dx*dy*dz` から導出すべき)。

### P5. メッシュの出自が記録されない

- `mesh_shifted_merged` / `mesh_refined_3x` がどの設定(細分化率・メッシュサイズ・
  ジオメトリ版)から生成されたか、リポジトリ内のデータとして残っていない。
- 実害の例: 旧メッシュは基板を y−1mm シフト、新ジオメトリは全体 y=+1mm。
  旧パルスUDFの `Y0=0.0` ハードコードは 3x メッシュでは**吸収体中心を外していた**。

### P6. GUI状態の混入

`role_map`(UUIDキー)、`elmer_overrides`、`autosave/` 群など、
シミュレーション定義ではない編集ツール由来の状態が主ファイルに同居している。

## 2. 設計原則

1. **宣言的な単一ソース**: 値は「式」1系統のみ。数値化はビルド時に行い、
   設定ファイルへの書き戻しは禁止(検証のみ)。
2. **このリポジトリで完結**: 単位付き式評価器(~150行)をベンダリングし、
   外部リポジトリ依存を切断。回路パラメータの取り込みは明示的なimportスクリプトで。
3. **生成と手書きを混在させない**: 実行されるSIFは全て `generated/cases/` に生成。
   手で編集するのは設定ファイルとテンプレートだけ。
4. **暗黙のデフォルト禁止**: UDFは必要な定数が無ければ `Fatal` で停止。
5. **メッシュはレジストリで管理**: 生成レシピと出自(PROVENANCE)をデータとして記録し、
   ケースはメッシュ名を参照する。派生量(離散正規化など)はビルド時に自動計算。

## 3. 目標アーキテクチャ

```
project.json            # 唯一の設定(下記スキーマ)
build.py                # 検証 → 派生量計算 → generated/cases/*.sif 生成
run.py <case>           # restartチェーン解決 → ElmerSolver → results/<case>/
build_mesh.py <mesh>    # gmsh → ElmerGrid → meshes/<name>/ + PROVENANCE.json
src/                    # UDF (.f90)
templates/              # SIFテンプレート (steady / transient / pulse / constant_power)
generated/              # ビルド生成物(全て再生成可能、手動編集禁止)
results/<case>/         # ソルバ出力(vtu/ep/result/CSV/ログ/manifest)
```

### スキーマ案(抜粋)

```jsonc
{
  "parameters": {                      // 式のみ。数値ブロックは廃止
    "abs_dx":  "1[mm]",
    "I_bias":  "715[uA]",
    "I_0":     "I_bias*R_sh/(R_0+R_sh)",
    "L_tes":   "12.3[nH]",             // ← 現在UDF内に埋没している値
    "TES_volume": "TES_Au_dx*TES_Au_dy*TES_dz"
  },
  "materials": { "Pb": { "rho": "9860", "cp": "3.26e-5", "k": "0.0168" }, ... },
  "meshes": {
    "shifted_merged": { "recipe": { "mesh_min": "50[um]", "refine": 1 }, "dir": "mesh_shifted_merged" },
    "refined_3x":     { "recipe": { "mesh_min": "50[um]", "refine": 3 }, "dir": "mesh_refined_3x" }
  },
  "cases": {
    "tes_steady_3x": {
      "template": "steady", "mesh": "refined_3x",
      "solver": { "nonlinear_max_iterations": 120, "tolerance": 1e-8 }
    },
    "tes_pulse_20ms_3x": {
      "template": "pulse", "mesh": "refined_3x",
      "restart_from": "tes_steady_3x",
      "pulse": { "energy": "1332[keV]", "start": "20.02[ms]",
                 "duration": "1[ns]", "sigma": "50[um]", "center": "auto" },
      "timesteps": [[ "1[ms]", 20 ], [ "18[us]", 1 ], [ "1[us]", 2 ], [ "1[ns]", 1 ],
                    [ "10[ns]", 10 ], [ "100[ns]", 9 ], [ "1[us]", 9 ], [ "10[us]", 9 ],
                    [ "100[us]", 9 ], [ "1[ms]", 79 ]],
      "solver": { "nonlinear_max_iterations": 25, "tolerance": 3e-7 }
    }
  }
}
```

### build.py が自動で行うこと

- スキーマ検証+式の評価(循環・未定義参照はエラー)。**書き戻さない**。
- 派生量: `I_0`、`TES_volume`、パルスの `center`(メッシュのabs体から重心を取得)、
  **`Pulse Discrete Norm`(メッシュ節点上のガウシアンをFE積分)** — 現在の手計算を自動化。
- テンプレート展開で `generated/cases/<case>.sif` を出力
  (Constantsブロック、タイムステップ配列、`Exec/Output Intervals` のステージ配列も生成)。
- `generated/parameter_summary.csv` の出力(現行機能の継続)。

### run.py が自動で行うこと

- `restart_from` の `.result` が無ければ先に依存ケースを実行(チェーン解決)。
- `Results Directory` を用いて出力を `results/<case>/` に隔離
  (メッシュディレクトリに出力が堆積する現状を解消。フェーズ2冒頭で挙動検証)。
- 実行ログ・終了コード・入力ハッシュ(project.json / メッシュ / DLL)を
  `results/<case>/manifest.json` に記録 → `freeze_repro_run.py` はこれを流用して簡素化。

### UDF の変更

- 全定数を必須化(`Found` チェック → 欠落時 `CALL Fatal`)。既定値の削除。
- `TES Inductance` を Constants 経由に(ハードコード廃止)。

## 4. 移行計画

| フェーズ | 内容 | 目安 | 解消する問題 |
|---|---|---|---|
| 0 | 式評価器のベンダリング(外部依存切断)/ reconcile を検証専用化(書き戻し・autosave廃止)/ `TES Inductance`・`TES_volume` の定数化 | 半日 | P1, P2(書き戻し), P4の一部 |
| 1 | ケースのデータ駆動化: 全8ケース(3x系含む)を `cases` セクション+テンプレートへ移行、discrete norm 自動計算、UDF定数必須化 | 1日 | P3, P4, P2(残り) |
| 2 | メッシュレジストリ+`build_mesh.py`+PROVENANCE、`run.py`(restartチェーン+`Results Directory`)、freeze統合、`role_map`等のGUI状態を別ファイルへ分離 | 1日 | P5, P6 |

### 検証(各フェーズ共通)

- ゴールデン回帰: `tes_steady_3x` → `tes_pulse_20ms_3x` を再実行し、
  `artifacts/series/tes_pulse_20ms_3x_series.csv`(現行確定結果)と全列を相対誤差
  1e-6 以内で一致確認。ベースメッシュ系は `case_tes_shunt_internal` の定常値で確認。
- フェーズ0完了時点で「外部リポジトリを一時リネームしてもビルドが通る」ことを確認。

## 5. 論点・リスク

- **設定ファイル形式**: コメント可能な TOML への移行が理想だが、
  Thermal-and-Electoric-Sim 側のGUI/ツールが同形式のJSONを読む場合は互換を壊す。
  → 既定案は「JSON継続+スキーマ整理」。GUI連携が不要なら Phase 1 で TOML 化を検討。
- **回路パラメータの同期**: `tes_test2.json`(外部リポジトリ)との整合は、
  自動同期ではなく `python scripts/prep/import_circuit_params.py`(明示実行・差分表示)
  とし、暗黙の依存を作らない。
- **`Results Directory` の挙動**(相対パス・ディレクトリ自動作成の有無)は
  Elmer 26.1 で要確認。問題があれば従来どおりメッシュディレクトリ出力+run.py が
  結果を `results/<case>/` へ移動する方式にフォールバック。
- 手書きSIFの廃止により、一時的な実験は「ケース定義の複製+パラメータ変更」で行う
  文化に変わる(実験の使い捨てSIFを作りたい場合は `generated/` 外で明示的に)。
