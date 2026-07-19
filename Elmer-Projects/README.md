# Elmer-Projects

このリポジトリは、TES の Elmer 熱計算を「編集する場所」と「生成・保存されるもの」を分けて扱う前提で整理しています。

## まず触る場所

- `elmer_project.json`: 現在の計算条件・形状定義を含む唯一の主ファイル
- `generate_project_geometry.py`: JSON から Gmsh/Elmer 用ジオメトリを生成
- `sync_elmer_parameters.py`: JSON から `generated/` の SIF 断片を再生成
- `tes_transient_heat_source.f90`: TES 電熱結合 UDF(定常/時間依存/パルス熱源)
- `tes_heat_source.f90`: 定電力熱源 UDF(constant power ケース用)

## ケース定義と実行

**ケースは `elmer_project.json` の `cases` セクションで定義**し、
`python sync_elmer_parameters.py` が自己完結な SIF を `generated/cases/` に生成します。
手書き SIF は存在しません(生成物を直接編集しない)。

定義済みケース(`template` / メッシュ):

- `case_tes_shunt_internal`(steady / base): 定常平衡 → `.result` 保存(pulse_1332kev の restart 元)
- `case_constant_power`, `case_constant_power_3x_refined`(steady): レガシー検証ケース
- `case_tes_steady_3x_refined`(steady / 3x): 定常平衡場 → `.result` 保存(pulse_20ms の restart 元)
- `case_tes_shunt_transient`, `case_tes_shunt_transient_3x_refined`(transient): 100 ms 時間依存
- `case_tes_pulse_1332kev`(pulse / base): t=0 から 0.5 µs 窓で 1332 keV
- `case_tes_pulse_20ms_3x_refined`(pulse / 3x): 20.02 ms に 1 ns 矩形窓で 1332 keV

パルス系はビルド時に**パルス中心(absの体積重心)と離散正規化係数をメッシュから自動計算**します。
値はケース定義に式(例 `"20.02[ms]"`, `"1332[keV]"`)で書けます。

### 2TES・長尺吸収体構成(dual-TES)

単ピクセル(absがTES基板上に載る)に加え、**1本の細長いabsの両端にTES基板スタックを
置く2TES構成**があります。電気回路はL/R独立の2系統で、UDFが3インスタンス化されています
(`TESTransientHeatSource` / `...L` / `...R`、定数プレフィックス無 / `TES L ` / `TES R `)。
形状は `geometries` レジストリで管理(`single_pixel` / `dual_tes`)し、各メッシュが
`geometry` を参照します。設計と実装記録は `docs/dual_tes_plan.md` を参照。

- `case_dual_steady`(steady / `mesh_dual_base`): 定常平衡 → `.result`(dual系のrestart元)
- `case_dual_shunt_transient`(transient): 無パルス過渡(ベースライン)
- `case_dual_pulse_center`(pulse): 中央入射(x≈0)。L/R対称応答
- `case_dual_pulse_offset`(pulse): 偏心入射(x=+1.5mm)。位置感応(R側が先行応答)

パルス中心はケース定義で `"center": "auto"`(abs重心)または
`{"x": "1.5[mm]", "y": "auto", "z": "auto"}` 形式の明示指定が可能です。
検証は `scripts/analysis/dual_series_analysis.py`(系列CSV)と
`dual_vtu_verify.py`(VTU場・エネルギー収支)。

実行は `run.py` 経由が標準です(restart 依存を自動解決し、出力を `results/<case>/` に
整理して `manifest.json` を書きます):

```powershell
elmerf90 tes_transient_heat_source.f90 -o tes_transient_heat_source_t0.dll
python run.py case_tes_pulse_20ms_3x_refined
```

依存ケース(`case_tes_steady_3x_refined`)の `.result` が無ければ先に自動実行されます。
`--dry-run` で実行計画のみ表示、`--force-deps` で依存も再実行。
`ElmerSolver generated\cases\<case>.sif` の直接実行も従来どおり可能です
(その場合の出力はメッシュディレクトリに残ります)。

メッシュは `meshes` レジストリで管理し、`python build_mesh.py <mesh名>` で
再生成できます(ジオメトリ生成器 gmsh 一式は `scripts/support/vendored/geometry/` に
ベンダリング済みで、外部リポジトリ依存はありません)。
各メッシュディレクトリの `PROVENANCE.json` にレシピとハッシュを記録します。

## ディレクトリ構成

- `docs/`: 説明資料(モデル定義・数値スキームの詳細は `docs/README_TES_elmer.md`)
- `generated/`: JSON から生成された SIF 断片(手動編集しない)
- `mesh_shifted_merged/`, `mesh_refined_3x/`: 単ピクセルメッシュ(ソルバ出力もここに書かれる)
- `mesh_dual_base/`: 2TES・長尺吸収体メッシュ(dual-TES 構成、19 body)
- `gmsh/`: ジオメトリ・メッシュ生成の中間物
- `scripts/analysis|visualization|prep|support/`: 抽出・可視化・準備・補助
- `artifacts/series/`: 確定した時系列 CSV(git 管理)/ `artifacts/plots/`: 図
- `reference/`: 外部由来の参照データ(COMSOL 時系列 `SignglePixel.txt`、`tes.json`, `tes_test2.json`)
- `runs/`: 凍結した再現用ラン(`python freeze_repro_run.py <name>`)
- `archive/`: 過去の実験ケース・結果(例: `cases_mortar_debug_202607/`)
- `legacy/`: 旧ケース・旧メッシュ(`legacy/meshes/` に引退メッシュ)
- `scripts/support/vendored/`: 単位付き式評価器(Thermal-and-Electoric-Sim からベンダリング)

## 基本ワークフロー

1. `elmer_project.json` を編集(**すべて式で記述**: `parameter_expressions`、materials の
   expression、geometry の `*_expr`。数値はビルド時に式から導出され、ファイルには持たない)
2. `python sync_elmer_parameters.py`(条件だけ変えた場合)/ `python build_mesh.py <mesh名>`(形状も変えた場合)
3. 必要なら `ElmerGrid ...` と `ElmerSolver ...` を実行
4. 結果確認は `scripts/analysis/` を使用
5. 確定結果は `python freeze_repro_run.py <run_name>` で `runs/` に凍結

## 代表コマンド

```powershell
python generate_project_geometry.py
python sync_elmer_parameters.py
ElmerSolver generated\cases\case_constant_power.sif
python scripts\analysis\summarize_tes_temperature.py
python freeze_repro_run.py current_reference
```

詳細は `docs/README_TES_elmer.md` を参照してください。
