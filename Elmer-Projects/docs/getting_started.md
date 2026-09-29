# シミュレーションの始め方

新しい計算は **１実験につき１つのケースファイル**を `cases/` に置きます。まずケースを選びます。

```powershell
python main.py list
python main.py show single_pixel_alpha240
python main.py run single_pixel_alpha240 --dry-run
python main.py run single_pixel_alpha240
```

`show` は使用するモデル、TES の数・位置、`alpha`、バイアスなどを表示します。`--dry-run` は生成・計算しません。`run` は必要なメッシュと回路ライブラリを用意し、定常解とパルス応答を計算します。

## alpha を変えた単ピクセル波形

[`cases/single_pixel.toml`](../cases/single_pixel.toml) をコピーして別名を付け、`[circuit]` に `alpha = 240` のように書きます。[`cases/single_pixel_alpha240.toml`](../cases/single_pixel_alpha240.toml) がその例です。`model = "single_pixel"` は変えません。元の h8 メッシュを再利用し、alpha ごとに別の定常解から波形を計算します。

```powershell
python main.py show single_pixel_alpha240
python main.py run single_pixel_alpha240
```

## 形状を変えるとき

ケースの `model` で既存のモデルを選びます。`single_pixel` は１ TES、`post_four_tes` は２ TES の部品を２組置く PoST 型です。既存の PoST を使うだけなら[そのケース](../cases/post_four_tes_nominal.toml)をコピーし、モデルファイルは編集しません。まだない配置・メッシュを作るときは `models/` に新しいモデルファイルを作り、ケースの `model` をその名前にします。

| 編集対象 | 内容 |
| --- | --- |
| `cases/<名前>.toml` | モデルの選択、`alpha`、バイアス、個別 TES の回路値、パルス、時間窓 |
| `models/<名前>.toml` | 共通値、材料、TES 部品、配置、メッシュ |

メッシュは形状を計算用の要素に分割したものです。物理値だけの変更では同じメッシュを使い、形状・要素寸法の変更では新しいメッシュを生成します。単ピクセルモデルは既存の `mesh_hybrid_fullconf_h8` を参照するため、その形状を変えるには別の生成モデルを作ります。

生成した内部 JSON と SIF は `generated/`、メッシュは `work/meshes/`、実行結果は `results/` に保存します。旧ケースは `python main.py list --all` で表示できます。実行には Gmsh、ElmerGrid、ElmerSolver、`elmerf90` が必要です。
