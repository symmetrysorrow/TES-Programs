# cases: legacy TOML

このディレクトリは `tes_sim.py` 系の旧 TOML を既存結果の再現用に残す場所です。
新しい計算は `projects/<名前>.toml` と対応する model TOML を作成し、
`python main.py run <名前>` で ElmerGPU Phase24 backend を使います。

```powershell
python main.py legacy2 list
python main.py legacy2 toml show projects/legacy/cases/<名前>.toml
```

旧形式は `legacy2` と従来の `legacy` コマンドから実行できます。
