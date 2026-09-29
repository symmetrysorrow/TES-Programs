# projects: 新しい計算の入口

1実験につき1つの `projects/<名前>.toml` を作ります。形状・材料・配置・mesh は
`projects/models/` の model TOML にまとめ、project から相対パスで参照します。

```toml
schema_version = 3
model = "models/single_pixel.toml"
```

回路値は override ではなく、TES ID ごとの最終値を project に書きます。

```toml
[circuits."pixel.tes"]
I_bias = "715[uA]"
R_sh = "3.9[mohm]"
L_tes = "12.3[nH]"
R_0 = "15.527[mohm]"
alpha = 256.46
beta = 5.03
```

TES ID は model の layout から決まり、`python main.py show <project>` で位置と一緒に
確認できます。

`templates/pipeline_template.json` は CLI が内部用 solver 設定を組み立てるための雛形です。
個別の計算条件は TOML に記述し、CLI が作る中間 JSON は `workspace/generated/` に保存します。
