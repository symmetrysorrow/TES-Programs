# Elmer-Projects のディレクトリとデータ

`Elmer-Projects/` を作業ディレクトリとした案内です。新しい計算は `projects/` から始めます。

| 段階 | 場所 | 中身・用途 |
| --- | --- | --- |
| 入力 | `projects/*.toml` (`schema_version = 3`) | 1実験の入口。model の相対パスと最終回路値、pulse、時間条件を持つ。 |
| 入力 | `projects/models/*.toml` | 再利用する形状、材料、TES配置、メッシュ。通常の条件変更では触らない。 |
| 旧方式の入力 | `cases/*.toml` | `tes_sim.py` 系の旧 TOML。既存結果の再現用。 |
| 旧方式の入力 | `elmer_project.json`、`elmer_project_*.json` | `run.py` 系の旧 JSON。 |
| 入力 | `config/` | 探索・解析などに使う設定ファイル。 |
| 入力資料 | `reference/` | COMSOL 時系列など比較元となる外部データ。 |
| 生成 | `gmsh/` | Gmsh のジオメトリ・メッシュ生成に伴うファイル。 |
| 生成 | `generated/`、`generated/cases/` | 実行用の内部 JSON / SIF。原則として手で編集しない。 |
| 作業中 | `work/meshes/` | 再生成可能な Elmer mesh と restart。 |
| 作業中 | `results/<case>/` | 実行出力。内部 solver case 名で保存される場合がある。 |
| 作業中 | `results/raw/`、`logs/` | 直接実行や調査の未整理出力。 |
| 確定 | `artifacts/<調査名>/` | 比較表、図、短い報告など。 |
| 実行記録・再現用 | `runs/<名前>/` | 旧 workflow の state / snapshot。 |
| 履歴 | `archive/` | 過去の実験・結果。 |
| 旧方式 | `legacy/` | 旧 SIF、旧 mesh、当時のログ。 |
| コード・説明 | `scripts/`、`tests/`、`docs/` | 実行・解析コード、検証、説明資料。 |

## project と model

通常ユーザーが最初に開くのは `projects/<名前>.toml` だけです。そのファイルに `model = "models/single_pixel.toml"` と書かれているため、形状や mesh の詳細が必要になったときだけ同じ `projects/` 配下の model を開きます。model を geometry / materials / mesh の別ファイルへさらに分割しません。

## project と結果

- project: 人が編集する計算条件
- model: 再利用する物理構造
- generated: project + model から作る内部入力
- work/results: 計算途中と結果
- artifacts: 残す価値のある解析結果

旧 `cases/` は新規計算の入口ではありません。新projectは `python main.py list`、旧形式は `python main.py legacy list` と入口を分けています。通常コマンドが旧形式へ自動フォールバックすることはありません。

## .dat / .csv

拡張子ではなく出所と用途で保存先を判断します。時系列や反復ログは対応する実行結果、確定した比較表は `artifacts/`、外部比較元は `reference/` に置きます。
