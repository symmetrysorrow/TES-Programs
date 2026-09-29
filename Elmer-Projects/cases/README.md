# cases: legacy TOML

このディレクトリは `tes_sim.py` 系の旧 TOML を既存結果の再現用に残す場所です。
新しい計算は `../projects/<名前>.toml` を作成し、`python main.py run <名前>` を使います。

```powershell
python main.py legacy list
python main.py legacy toml show cases/<名前>.toml
```

通常の `python main.py list/show/run` は `projects/*.toml` 専用です。このディレクトリの入力は `legacy` 名前空間からだけ案内します。
