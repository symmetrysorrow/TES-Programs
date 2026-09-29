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
python main.py show single_pixel_alpha240
python main.py run single_pixel_alpha240 --dry-run
python main.py run single_pixel_alpha240
```

`show` は参照する model、TES の ID と位置、各 TES の最終的な `I_bias`、`R_sh`、`alpha` などを表示します。`--dry-run` は生成・計算を行いません。

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

`projects/single_pixel.toml` をコピーして別名を付け、その TES の `[circuits."pixel.tes"]` を編集します。

```toml
[circuits."pixel.tes"]
I_bias = "715[uA]"
R_sh = "3.9[mohm]"
L_tes = "12.3[nH]"
R_0 = "15.527[mohm]"
alpha = 240
beta = 5.03
```

これは override ではなく、この計算で使う**最終値**です。model、component、group から同名の回路値を継承して上書きする仕組みはありません。複数 TES の project では、model が定義する TES ID ごとに1つずつ完全な `circuits` テーブルを書きます。

ID が分からない場合は `python main.py show <名前>` で位置と一緒に確認できます。

## 形状・材料・メッシュを変える

形状・材料・メッシュは model にまとまっています。通常の bias / alpha / pulse 変更で model を編集する必要はありません。新しい形状を作る場合は `projects/models/` に新しい model TOML を作り、project の `model` パスをそのファイルへ向けます。

単ピクセル model は検証済みの `mesh_hybrid_fullconf_h8` を参照しています。この prebuilt mesh と矛盾する寸法変更はエラーになります。新しい寸法を使う場合は、生成型 mesh を持つ別 model を作成します。

生成した内部 JSON と SIF は `workspace/generated/`、メッシュは `workspace/work/meshes/`、実行結果は `outputs/results/` に保存します。旧 TOML / JSON を再現する場合だけ `python main.py legacy list` を使います。通常の `list/show/check/mesh/run` は新project専用です。
