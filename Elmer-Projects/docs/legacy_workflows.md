# Legacy workflows

この文書は、既存結果の再現のために残している旧 JSON / 旧 TOML 経路の説明です。
**新しい計算は `projects/*.toml` と通常の `main.py list/show/check/mesh/run` を使ってください。**

旧経路は通常CLIと混ざらないよう、`legacy` 名前空間からだけ案内します。

```powershell
python main.py legacy list
python main.py legacy check
python main.py legacy show <case-or-mesh>
python main.py legacy mesh <mesh>
python main.py legacy run <case> --dry-run
python main.py legacy toml show projects/legacy/cases/<name>.toml
```

旧 JSON を別ファイルから読む場合だけ、legacyサブコマンドの `--project` を使います。

## 旧 JSON 経路のケース定義と実行

**`run.py` で使うケースは `projects/legacy/elmer_project.json` の `cases` セクションで定義**し、
`python src/sync_elmer_parameters.py` が自己完結な SIF を `workspace/generated/cases/` に生成します。
手書き SIF は存在しません(生成物を直接編集しない)。
別系統の `tes_sim.py` は `projects/legacy/cases/*.toml` を入力として使います。

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
検証は `src/analysis/dual_series_analysis.py`(系列CSV)と
`dual_vtu_verify.py`(VTU場・エネルギー収支)。

実行は `run.py` 経由が標準です(restart 依存を自動解決し、出力を `outputs/results/<case>/` に
整理して `manifest.json` を書きます):

```powershell
elmerf90 src/fortran/tes_transient_heat_source.f90 -o build/udf/tes_transient_heat_source_t0.dll
python src/run.py case_tes_pulse_20ms_3x_refined
```

依存ケース(`case_tes_steady_3x_refined`)の `.result` が無ければ先に自動実行されます。
`--dry-run` で実行計画のみ表示、`--force-deps` で依存も再実行。
`ElmerSolver workspace/generated\cases\<case>.sif` の直接実行も従来どおり可能です
(その場合の出力はメッシュディレクトリに残ります)。

### GPU / AMGX

GPU起動スクリプトとAMGX設定は削除しました。ここに記載していたGPU実行手順は現在利用できません。

メッシュは `meshes` レジストリで管理し、`python src/build_mesh.py <mesh名>` で
再生成できます(ジオメトリ生成器 gmsh 一式は `src/support/vendored/geometry/` に
ベンダリング済みで、外部リポジトリ依存はありません)。
各 `workspace/work/meshes/<mesh名>/PROVENANCE.json` にレシピとハッシュを記録します。

## 旧 JSON 経路の基本ワークフロー

1. `projects/legacy/elmer_project.json` を編集(**すべて式で記述**: `parameter_expressions`、materials の
   expression、geometry の `*_expr`。数値はビルド時に式から導出され、ファイルには持たない)
2. `python src/sync_elmer_parameters.py`(条件だけ変えた場合)/ `python src/build_mesh.py <mesh名>`(形状も変えた場合)
3. 必要なら `ElmerGrid ...` と `ElmerSolver ...` を実行
4. 結果確認は `src/analysis/` を使用
5. 比較・診断の確定結果は `outputs/analysis/<調査名>/` に記録。旧 `case_constant_power` の再現用一式だけは、必要な旧ファイルがある場合に `freeze_repro_run.py` で `outputs/runs/` に凍結可能

## 代表コマンド

```powershell
python src/generate_project_geometry.py
python src/sync_elmer_parameters.py
ElmerSolver workspace/generated\cases\case_constant_power.sif
python src\analysis\summarize_tes_temperature.py
```

詳細は `docs/README_TES_elmer.md` を参照してください。
