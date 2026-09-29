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
python main.py show single_pixel_alpha240
python main.py run single_pixel_alpha240 --dry-run
python main.py run single_pixel_alpha240
```

新しい計算の入口は **`projects/<実験名>.toml`** です。project は `model = "models/single_pixel.toml"` のように、同じ `projects/` 配下の model を明示的な相対パスで参照します。

- `projects/*.toml`: その計算で実際に使う回路値、パルス、時間窓、solver 条件
- `projects/models/*.toml`: 形状、材料、TES 配置、メッシュ
**project/model 間の汎用 override はありません。** `I_bias`、`R_sh`、`alpha` などの実行時回路値は model に置かず、project の `[circuits."<TES ID>"]` に各 TES の最終値を明示します。ID は安定した参照名としてだけ使い、component → group → project → TES のような値の上書きチェーンは作りません。

メッシュ条件は複雑で通常は触らないため model 側にまとめています。project を読めば「どの model を使うか」は分かり、詳細が必要なときだけその1ファイルを開けば済みます。`python main.py show <project>` で model パス、TES ID・位置、実効回路値を確認できます。

詳しい編集手順は [シミュレーションの始め方](docs/getting_started.md)、構成方針は [project/model 設計](docs/toml_unification_plan.md) を参照してください。

## まず触る場所

- `projects/*.toml`: 1実験ごとの最終計算条件
- `projects/models/*.toml`: 再利用する形状・材料・配置・メッシュ
- `main.py`: `list / show / check / mesh / run` の通常入口
- `docs/getting_started.md`: 新しい project の編集手順

旧 JSON / `projects/legacy/cases/*.toml` を再現する必要がある場合だけ、
[Legacy workflows](docs/legacy_workflows.md) を参照してください。通常CLIから旧形式への自動フォールバックはありません。

## Legacy

旧方式は明示的に `legacy` 名前空間へ分離しています。

```powershell
python main.py legacy list
python main.py legacy show <旧JSONのcaseまたはmesh>
python main.py legacy run <旧JSONのcase>
python main.py legacy toml show projects/legacy/cases/<名前>.toml
```

新規計算でこれらのコマンドを使う必要はありません。詳細は
[Legacy workflows](docs/legacy_workflows.md) にまとめています。

## ディレクトリ構成

入力から結果までの置き場所、`archive/` と `legacy/` の区別、`.dat` / `.csv` の読み方は
[ディレクトリとデータの案内](docs/directory_guide.md) にまとめています。
モデル定義・数値スキームの詳細は [Elmer の技術資料](docs/README_TES_elmer.md) を参照してください。
