# project / model TOML 設計

日付: 2026-09-29 / 状態: 実装済み

## 方針

新しいユーザー向け入口を `projects/*.toml` に統一し、その配下の `projects/models/*.toml` に再利用する物理モデルを置く。

```text
projects/
├── single_pixel.toml
├── single_pixel_alpha240.toml
├── post_four_tes_nominal.toml
└── models/
    ├── single_pixel.toml
    └── post_four_tes.toml
```

project は `model = "models/single_pixel.toml"` のように明示的な相対パスで model を参照する。これにより `projects/legacy/cases/` と `models/` がルートで並列に見える問題をなくし、「project が入口、その下に model がある」という位置関係をファイル構造にも反映する。

## 所有権

model が所有するもの: geometry / 寸法、materials、component と TES ID、layout / 配置、mesh。

project が所有するもの: 各 TES の最終回路値、pulse、time window / timesteps、solver の実行条件、field 保存の有無。

同じ実行時回路キーを model と project の両方に置かない。model に `I_bias`, `R_sh`, `L_tes`, `R_0`, `alpha`, `beta`, `I_0` が存在した場合は入力エラーにする。

## override を廃止

旧新形式では model parameter → component → group → case → individual TES の順に回路値を上書きできた。この多段 override は廃止する。

project は model が生成する TES ID ごとに `[circuits."<TES ID>"]` を持ち、`I_bias`, `R_sh`, `L_tes`, `R_0`, `alpha`, `beta` の最終値をすべて明示する。不足、余分なキー、存在しない TES ID は実行前にエラーにする。`I_0` はこれらから内部生成時に導出する。

ID は値の継承階層ではなく、TES個体を安定して参照するためだけに使う。

## 互換性

旧 `projects/legacy/cases/*.toml` / `tes_sim.py` と `elmer_project.json` / `run.py` は再現用として残す。新式 project/model TOML の通常 `run` は Elmer WSL backend を使い、以前のローカル実行は `main.py legacy1 ...`、旧 JSON / 旧ケース TOML は `main.py legacy2 ...` と従来の `main.py legacy ...` に残す。通常コマンドは旧入力形式へフォールバックしない。

## Public TES IDs

Model format v3 maps each component-local TES ID to a stable public ID in `layout.groups[].tes_ids`. Projects use public IDs as `[circuits.tes_1]` keys. Model format v2 remains supported with its historical `<group>.<local_id>` IDs. Mesh identity depends on TES positions and physical settings, not public ID spelling.
