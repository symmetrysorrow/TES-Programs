# Elmer-Projects

Source and compiled library locations are described in [docs/source_layout.md](docs/source_layout.md).

このリポジトリは、TES の Elmer 熱計算を「編集する場所」と「生成・保存されるもの」を分けて扱う前提で整理しています。

**どこに何があるか迷ったら [ディレクトリとデータの案内](docs/directory_guide.md) を参照してください。**
計算データは `workspace/`（作業データ）と `outputs/`（計算結果・再現用実行）にまとめています。旧形式は [Legacy workflows](docs/legacy_workflows.md) に分離しています。
`tools/` はローカル実行用の Elmer 配布物です。別途 Elmer をインストールしている場合は使わずに実行できます。

## ここから始める

`Elmer-Projects/` で次を実行します。

```powershell
python -m pip install questionary
python main.py
```

メニューから操作と project を選べます。実行時は dry-run（計画の確認）か計算開始かを選びます。従来のコマンド形式も引き続き使えます。

```powershell
python main.py list
python main.py run single_pixel_alpha240 --dry-run
python main.py run single_pixel_alpha240
```

`main.py run <project>` keeps the project/model TOML interface and executes the
steady solve and transient with the ElmerGPU Phase24 backend in WSL. The
validated fast solver settings are applied to the generated internal case;
the TOML remains the source for geometry, circuit, pulse, and time-window
values.

Use `python main.py legacy1 run <project>` to reproduce the former local Elmer
route. The older JSON and `src/tes_sim.py` workflows remain under
`python main.py legacy2 ...` (`python main.py legacy ...` is kept as an alias).

新しい計算の入口は **`projects/<実験名>.toml`** です。project は `model = "models/single_pixel.toml"` のように、同じ `projects/` 配下の model を明示的な相対パスで参照します。

- `projects/*.toml`: その計算で実際に使う回路値、パルス、時間窓、solver 条件
- `projects/models/*.toml`: 形状、材料、TES 配置、メッシュ
**project/model 間の汎用 override はありません。** `I_bias`、`R_sh`、`alpha` などの実行時回路値は model に置かず、project の `[circuits."<TES ID>"]` に各 TES の最終値を明示します。ID は安定した参照名としてだけ使い、component → group → project → TES のような値の上書きチェーンは作りません。


詳しい編集手順は [シミュレーションの始め方](docs/getting_started.md)、構成方針は [project/model 設計](docs/toml_unification_plan.md) を参照してください。

## まず触る場所

- `projects/*.toml`: 1実験ごとの最終計算条件
- `projects/models/*.toml`: 再利用する形状・材料・配置・メッシュ
- `docs/getting_started.md`: 新しい project の編集手順

以前の計算経路や旧入力を再現する場合は、
[Legacy workflows](docs/legacy_workflows.md) を参照してください。通常CLIから旧入力形式への自動フォールバックはありません。

## Legacy

旧方式は互換用のコマンドに残しています。

```powershell
python main.py legacy1 run single_pixel_alpha240 --dry-run
python main.py legacy2 list
python main.py legacy2 toml show projects/legacy/cases/<名前>.toml
```

`legacy1` は従来の project/model TOML のローカル Elmer 実行、`legacy2` は旧 JSON と
`src/tes_sim.py` の入力を使う経路です。`legacy` コマンドも旧 JSON / TOML の互換入口として残ります。詳細は
[Legacy workflows](docs/legacy_workflows.md) にまとめています。

## ディレクトリ構成

入力から結果までの置き場所、`archive/` と `legacy/` の区別、`.dat` / `.csv` の読み方は
[ディレクトリとデータの案内](docs/directory_guide.md) にまとめています。
モデル定義・数値スキームの詳細は [Elmer の技術資料](docs/README_TES_elmer.md) を参照してください。

## Public TES IDs

Models assign stable project-facing TES IDs through `layout.groups[].tes_ids`; project circuits use keys such as `[circuits.tes_1]`. Edit the project TOML for run and circuit values, and the referenced model TOML for geometry, materials, layout mappings, and mesh settings.
