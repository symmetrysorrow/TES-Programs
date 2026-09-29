# runs: 実行状態と固定スナップショット

ここには二種類のデータがあります。

| 作成元 | 置かれるもの | 意味 |
| --- | --- | --- |
| `tes_sim.py` | `runs/<ケース名>/state.json`、`summary.json` | 実行状態と要約。入力・メッシュ・結果のコピーではない。 |
| `freeze_repro_run.py` | `runs/<名前>/manifest.json`、`README.md`、入力・生成物・メッシュ・結果、対応する ZIP | 旧 `case_constant_power` の再現用に固定した一式。 |

`freeze_repro_run.py` は旧定電力ケース専用です。ルートの `elmer_project.json`、`case_constant_power.sif`、旧メッシュ `mesh_shifted_merged/` などがそろっている環境でのみ使います。任意の `run.py` / `tes_sim.py` ケースを凍結する汎用機能ではありません。

```powershell
python freeze_repro_run.py current_reference
```

通常の実行結果は `results/<case>/`、軽量な比較結果は `artifacts/`、手書きの TOML 入力は `cases/` を参照してください。全体の使い分けは [ディレクトリとデータの案内](../docs/directory_guide.md) にあります。
