# cases: 計算条件

新しい計算では１実験につき１ファイルを作ります。`schema_version = 2` と `model = "single_pixel"` または `model = "post_four_tes"` を指定し、回路・パルスなどの条件を書きます。例は `single_pixel_alpha240.toml` です。`python main.py show single_pixel_alpha240` で確認し、`python main.py run single_pixel_alpha240` で実行します。モデルの中身は `../models/` にあります。

`single_pixel_h8.toml` などの旧形式は既存結果の再現用に残しています。`python main.py run cases/single_pixel_h8.toml` で実行できます。

`run.py` 系のケース定義はルートの `elmer_project.json` にあり、生成された SIF は `generated/cases/` に入ります。実行結果は `results/<case>/` に入ります。

全体の置き場所は [ディレクトリとデータの案内](../docs/directory_guide.md) を参照してください。
