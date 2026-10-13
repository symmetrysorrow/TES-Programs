# Elmer-Projects

TES の Elmer 計算プロジェクトです。日常の計算では `projects/*.toml` の条件と
`projects/models/*.toml` の形状・材料・メッシュ設定を編集し、`ElmerGPU` という名前の WSL
ディストリビューション上でソルバーを実行します。

ソースとコンパイル済みライブラリの配置は
[source layout](docs/source_layout.md)、生成物・結果の配置は
[directory guide](docs/directory_guide.md) を参照してください。

## CPU専用の導入（外部研究者向け・試験段階）

Windows 11 + Ubuntu 24.04 WSL2 がある環境では、GPUを使わずに
改造版 Elmer/HYPRE と必要なPythonライブラリを入れるための初回セットアップを
用意しました。詳細と制限事項は [CPU setup guide](docs/cpu_onboarding.md) を参照してください。

```powershell
python main.py setup
.\.venv\Scripts\python.exe main.py run single_pixel_generated --dry-run
.\.venv\Scripts\python.exe main.py run single_pixel_generated
```

`single_pixel_generated` はTOMLからメッシュを作る**未検証の試験用**設定です。
検証済みの `single_pixel_alpha240`（既存 h8 メッシュ）は変更していません。
生成メッシュとCOMSOL/既存h8の数値一致は、別途検証が必要です。
`setup` 自体も第三者環境での動作検証前です。

## 標準計算経路

```text
python main.py run <project>
  → scenario_cli.py (TOML 読み込み、ケース生成、restart 判定)
  → elmer_wsl_backend.py (設定適用、MPI 分割、WSL 起動)
  → WSL: src/run.py (SIF、manifest、結果・restart 管理)
  → ElmerSolver
```

`elmer_wsl` は Windows から WSL の Elmer 実行環境を選ぶバックエンドです。
これは HYPRE の GPU 演算を強制する指定ではありません。既定の
`iterative_hypre_pcg_boomeramg` は CPU HYPRE を使います。
`linear_system` を `_gpu` で終わる対応値にした場合に限り、生成 SIF の
`HYPRE GPU = Logical True` が有効になります。WSL 環境に GPU 対応
ElmerSolver/HYPRE がインストールされている必要があります。

旧 `--backend elmergpu` も互換エイリアスとして使えます。現在の WSL ディストリ名は既定で `ElmerGPU` です。新しいCPU専用環境は
`/opt/tes-elmer-cpu` に導入します（旧GPU環境は上書きしません）。

## 実行

対話メニューを使う場合は `questionary` をインストールして `python main.py`
を実行します。通常の実験は以下のコマンドでも操作できます。

```powershell
python -m pip install questionary
python main.py
python main.py list
python main.py run single_pixel_alpha240 --dry-run
python main.py run single_pixel_alpha240
python main.py run single_pixel_alpha240 --mpi-procs 4
```

`--dry-run` はメッシュ生成、MPI 分割、設定書き換え、ソルバー起動を行いません。
通常の実行では、標準計算プロジェクトに対して必要な MPI 分割を再利用または
ElmerGrid で生成してから計算します。

`projects/<name>.toml` は回路・パルス・ソルバー条件を、参照先モデル TOML は
形状・材料・TES 配置・メッシュ条件を管理します。詳しくは
[計算開始ガイド](docs/getting_started.md) を参照してください。
`project/model` 間に汎用の上書き連鎖はありません。TES ごとの回路値は project の
`[circuits."<TES ID>"]` に記述します。

## Legacy と検証資料

過去のローカル Elmer 実行と旧 JSON / TOML 入力は、再現用に
`python main.py legacy1 ...`、`python main.py legacy2 ...` として残しています。
これらは標準経路から独立しています。過去の実験ランナー、記録、Elmer の
ユーザー定義 SIF キーワードは互換性と再現性のため保持されています。
`Phase24` などの接頭辞を持つ一部キーワードは ElmerSolver 側との内部
インターフェースであり、Python の表示名とは別に維持しています。

ディレクトリ構成は [directory guide](docs/directory_guide.md)、旧手順は
[legacy workflows](docs/legacy_workflows.md) を参照してください。

モデルと数値スキームの補足は [Elmer TES technical notes](docs/README_TES_elmer.md)、
project/model の構成方針は [TOML design notes](docs/toml_unification_plan.md)
を参照してください。

## Public TES IDs

モデルは `layout.groups[].tes_ids` で安定した TES ID を割り当てます。
回路値は `[circuits.tes_1]` のようなキーで対応付けます。実行条件は project TOML、
形状・材料・配置・メッシュ条件は参照先の model TOML で編集してください。
