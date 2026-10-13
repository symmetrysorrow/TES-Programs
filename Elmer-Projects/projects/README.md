# projects: 新しい計算の入口

1実験につき1つの `projects/<名前>.toml` を作ります。形状・材料・配置・mesh は
`projects/models/` の model TOML にまとめ、project から相対パスで参照します。

```toml
schema_version = 3
model = "models/single_pixel.toml"
```

回路値は override ではなく、TES ID ごとの最終値を project に書きます。

```toml
[circuits.tes_1]
I_bias = "715[uA]"
R_sh = "3.9[mohm]"
L_tes = "12.3[nH]"
R_0 = "15.527[mohm]"
alpha = 256.46
beta = 5.03
```

確認できます。

`templates/pipeline_template.json` は CLI が内部用 solver 設定を組み立てるための雛形です。
個別の計算条件は TOML に記述し、CLI が作る中間 JSON は `workspace/generated/` に保存します。

通常の `python main.py run <project>` は内部 JSON を作ったあと、Elmer WSL backend
で計算します。以前のローカル Elmer 実行は `python main.py legacy1 run <project>` から使えます。

Public TES IDs are explicit model mappings. For example, `tes_ids = { tes = "tes_1" }` maps a component-local TES ID to the project-facing ID. Projects refer to `[circuits.tes_1]`. Edit the project TOML for circuit values and the model TOML for geometry and mesh settings.
