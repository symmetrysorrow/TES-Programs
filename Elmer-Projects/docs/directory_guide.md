# Elmer-Projects のディレクトリとデータ

`Elmer-Projects/` を作業ディレクトリとしたときの案内です。まず入力を選び、生成・実行してから、必要な結果だけを保存します。

| 段階 | 場所 | 中身・用途 |
| --- | --- | --- |
| 入力 | `cases/*.toml` (`schema_version = 2`) | 新しい計算の１実験ごとの条件。使うモデル名も指定する。 |
| 入力 | `models/*.toml` | 再利用する TES の形状、配置、メッシュ、共通物性。 |
| 旧方式の入力 | `elmer_project.json`、`elmer_project_*.json` | `run.py` 系の計算条件、形状、ケース定義。派生 JSON は特定の比較・実験用。 |
| 旧方式の入力 | `cases/*.toml` (バージョン指定なし) | `tes_sim.py` 系のケース入力。新形式とは `schema_version` で区別する。 |
| 入力 | `config/` | 探索・解析などに使う設定ファイル。 |
| 入力資料 | `reference/` | COMSOL の時系列など、比較元となる外部由来のデータ。実行結果の置き場ではない。 |
| 生成 | `gmsh/` | Gmsh のジオメトリ・メッシュ生成に伴うファイル。 |
| 生成 | `generated/`、`generated/cases/` | JSON から作る SIF 断片や実行用 SIF。`sync_elmer_parameters.py` で再生成するため、原則として手で編集しない。 |
| 作業中 | `work/meshes/` | 再生成できる Elmer メッシュと、そのメッシュに結び付く restart `.result`。`PROVENANCE.json` にレシピとハッシュを記録する。Git 管理外。 |
| 作業中 | `results/<case>/` | `run.py` / `tes_sim.py` のケース別実行出力。`run.py` は `manifest.json` も書く。Git 管理外。 |
| 作業中 | `results/raw/`、`logs/` | 直接実行や調査で得た未整理の出力・ログ。結果の意味を確認してから保存先を決める。 |
| 確定 | `artifacts/<調査名>/` | 比較表、図、短い報告など、後から参照したい軽量な解析結果。Git 管理対象。 |
| 実行記録・再現用 | `runs/<名前>/` | `tes_sim.py` の `state.json` / `summary.json`、または `freeze_repro_run.py` が作る固定スナップショット。用途は中のファイルで見分ける。Git 管理外。 |
| 履歴 | `archive/` | 日付や調査名でまとめて保存した過去の実験・結果。現在の入力として使わないもの。 |
| 旧方式 | `legacy/` | 旧 SIF、旧メッシュ、当時のログなど、以前の実装方式を調べるための資料。 |
| コード・説明 | `scripts/`、`tests/`、`docs/` | 実行・解析コード、検証コード、説明資料。 |

## `cases` と `runs` の違い

- `cases/<名前>.toml` (`schema_version = 2`): 新しいケース定義。例: `python main.py run single_pixel_alpha240`。
- `elmer_project.json` の `cases`: 旧 `run.py` 系のケース定義。`generated/cases/*.sif` はその生成物。
- 旧 `cases/*.toml`: `tes_sim.py` 系のケース入力。例: `python main.py run cases/single_pixel_h8.toml`。
- `results/<case>/`: 実行のたびに増える作業用出力。
- `runs/<名前>/state.json` / `summary.json`: `tes_sim.py` の実行状態と要約。計算一式のコピーではない。
- `runs/<名前>/manifest.json`: `freeze_repro_run.py` による過去の定電力ケースの固定スナップショット。

## `.dat` / `.csv` を見つけたら

拡張子だけでは保存場所を決められません。**出所と用途**で判断します。

| ファイル | 読み方・保存先 |
| --- | --- |
| `<case>_series.csv` | 時系列。ケースの実行出力として `results/<case>/` を見る。 |
| `<case>_iterations.csv` | ソルバー反復の診断値。同じケースの `results/<case>/` を見る。 |
| `*_a.dat`、`*_b.dat`、`*_sizes.dat` | 行列・右辺・サイズなどの数値診断用ダンプ。対応する調査の生データとして扱う。 |
| 比較表・集計 CSV | 結論と対応付けて残すものは `artifacts/<調査名>/` に置く。 |
| COMSOL など外部由来の値 | 比較元として使う入力は `reference/` に置き、出所を記録する。 |

ルート直下にも過去の直接実行による `.dat` / `.csv` が残っています。名前から内容を推測して一括移動すると、既存スクリプトの相対パスや手元の調査に影響します。新しい調査ではケース別の `results/` を使い、確定した集計だけを `artifacts/` に置きます。過去のファイルを整理するときは、参照しているコード・文書を確認し、調査単位で移します。

## 選び方の例

- 条件を変える → `cases/<名前>.toml` をコピーして別名のケースを作る。
- 既存の形状を切り替える → ケースの `model` を選び直す。新しい配置やメッシュを作る → `models/` に新モデルを作る。
- 計算結果を確認する → `results/<case>/` を開く。
- 比較結果を共有する → `artifacts/<調査名>/` に表と短い説明をまとめる。
- 旧定電力ケースの結果一式を固定する → 必要な旧ファイルがそろっている場合に `python freeze_repro_run.py <名前>` を使う。他のケースは実行条件と結果を調査単位で `artifacts/` に記録する。
- 古い実装を調べる → `legacy/`、過去の調査を調べる → `archive/` を開く。
