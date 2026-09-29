# Elmer-Projects のディレクトリとデータ

新しい計算は `projects/` から始めます。計算データはルート直下に散らばらないよう、`workspace/` と `outputs/` の2つにまとめています。

| 段階 | 場所 | 中身・用途 |
| --- | --- | --- |
| 入力 | `projects/*.toml` (`schema_version = 3`) | 1実験の入口。model の相対パスと最終回路値、pulse、時間条件を持つ。 |
| 入力 | `projects/models/*.toml` | 再利用する形状、材料、TES配置、メッシュ。 |
| 旧方式の入力 | `projects/legacy/cases/*.toml`、`elmer_project*.json` | 旧 workflow の再現用入力。 |
| 入力資料 | `docs/reference-data/` | COMSOL時系列など比較元となる外部データ。 |
| 作業データ | `workspace/work/` | Elmer mesh、restart、作業中の中間データ。 |
| 作業データ | `workspace/generated/` | project/model から生成する内部 JSON、SIF、図。 |
| 作業データ | `workspace/gmsh/` | Gmsh の形状・メッシュと生成途中のデータ。 |
| 計算出力 | `outputs/results/<case>/` | 計算ごとの結果、ログ、manifest。 |
| 再現用実行 | `outputs/runs/<name>/` | 固定した実行一式と再現用 ZIP。 |
| 確定資料 | `artifacts/<調査名>/` | 残す価値のある比較表、図、報告。 |
| 履歴 | `archive/`、`legacy/` | 過去の実験記録と旧 workflow。 |
| 実行環境 | `tools/` | ローカルの Elmer 配布物。Elmer を別途用意している場合は不要。Git 管理外。 |
| コード・説明 | `src/`、`tests/`、`docs/` | 実装、補助ツール、検証、説明資料。 |

## project と model

通常ユーザーが最初に開くのは `projects/<名前>.toml` です。そのファイルが同じ `projects/` 配下の model を参照します。形状や mesh の詳細を変更するときだけ model を開きます。

## workspace と outputs

- `workspace/` は再生成できる作業用データをまとめます。`work/` はElmer meshとrestart、`generated/` は内部入力、`gmsh/` は形状・メッシュ生成物です。
- `outputs/` は計算出力と再現用スナップショットをまとめます。`results/` はケースごとの通常出力、`runs/` は固定した再現パッケージです。
- `artifacts/` は比較・解析結果のうち、保存して共有するものを置きます。

旧 `projects/legacy/cases/` や JSON は新規計算の入口ではありません。新projectは `python main.py list`、旧形式は `python main.py legacy list` で扱います。

## .dat / .csv

拡張子ではなく出所と用途で保存先を判断します。時系列や反復ログは対応する `outputs/results/`、確定した比較表は `artifacts/`、外部比較元は `docs/reference-data/` に置きます。
