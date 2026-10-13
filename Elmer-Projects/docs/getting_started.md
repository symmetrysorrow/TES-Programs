# シミュレーションの始め方

新しい計算は **1実験につき1つの project TOML** を `projects/` に置きます。

コマンドを覚えずに始める場合は、`questionary` を入れてメニューを開きます。

```powershell
python -m pip install questionary
python main.py
```

メニューで project と操作を選べます。`計算を実行する` を選ぶと、dry-run で計画を確認するか、そのまま計算するかを選択できます。`python main.py --interactive` でもメニューを開けます。

スクリプトやラッパーから使う場合は、引き続き次のコマンド形式を使えます。

```powershell
python main.py list
python main.py run single_pixel_alpha240 --dry-run
python main.py run single_pixel_alpha240
```

`run` は Elmer WSL backend を使います。計算機には `ElmerGPU` WSL distribution と
CPU環境では `/opt/tes-elmer-cpu` の実行環境が必要です。以前のローカル Elmer 経路を再現する場合は
`python main.py legacy1 run single_pixel_alpha240` を使えます。旧 JSON / 旧ケース TOML
形式は `python main.py legacy2 ...` に残しています。


既定の backend 名は `elmer_wsl` です。旧 `elmergpu` 指定も互換エイリアスとして
受け付けます。`ElmerGPU` は現在の WSL ディストリの固有名であり、GPU 演算の指定では
ありません。CPU/GPU の選択は `[solver] linear_system` によって行います。

初回のCPU専用環境構築と自動メッシュ生成の試験例は
[CPU setup guide](cpu_onboarding.md) を参照してください。現行の
`single_pixel_alpha240` は既存の固定メッシュを必要とし、
`single_pixel_generated` は生成メッシュの検証前の試験用です。

## project と model

project は実行条件、model は通常変更しない構造条件です。

| 編集対象 | 内容 |
| --- | --- |
| `projects/<名前>.toml` | model の選択、各 TES の回路値、パルス、時間窓、solver |
| `projects/models/<名前>.toml` | 寸法、材料、TES 部品、配置、メッシュ |

project の `model` は名前ではなく相対パスです。

```toml
schema_version = 3
model = "models/single_pixel.toml"
```

これにより project を開いた時点で model の場所が分かります。

## 回路値を変える

`projects/single_pixel.toml` をコピーして別名を付け、その TES の `[circuits.tes_1]` を編集します。

```toml
[circuits.tes_1]
I_bias = "715[uA]"
R_sh = "3.9[mohm]"
L_tes = "12.3[nH]"
R_0 = "15.527[mohm]"
alpha = 240
beta = 5.03
```

これは override ではなく、この計算で使う**最終値**です。model、component、group から同名の回路値を継承して上書きする仕組みはありません。複数 TES の project では、model が定義する TES ID ごとに1つずつ完全な `circuits` テーブルを書きます。


## 形状・材料・メッシュを変える

形状・材料・メッシュは model にまとまっています。通常の bias / alpha / pulse 変更で model を編集する必要はありません。新しい形状を作る場合は `projects/models/` に新しい model TOML を作り、project の `model` パスをそのファイルへ向けます。

単ピクセル model は検証済みの `mesh_hybrid_fullconf_h8` を参照しています。この prebuilt mesh と矛盾する寸法変更はエラーになります。新しい寸法を使う場合は、生成型 mesh を持つ別 model を作成します。

## Public TES IDs

Use `[circuits.tes_1]` in a project. The referenced model assigns `tes_1` through `layout.groups[].tes_ids`; public IDs are independent of component-local and layout group IDs. Edit the project TOML for circuit and run values, and the model TOML for geometry, mesh settings, and ID mappings.

## Rebuilding the Elmer solver in WSL

The pinned Elmer/HYPRE source commits, source patches, Ubuntu package versions,
and build instructions are documented in [`elmergpu_build.md`](elmergpu_build.md).
