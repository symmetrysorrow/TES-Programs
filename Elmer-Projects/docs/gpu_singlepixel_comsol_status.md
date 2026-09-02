# GPU SinglePixel 比較状況

2026-09-02 時点で、RTX 3060 Ti の CUDA/AMGX 実行基盤と、COMSOL のパルス立上がりを解像する比較ケースを用意した。

## 再現ケース

```powershell
python scripts/prep/prepare_gpu_hybrid_prism_phase19.py
.scriptsun_singlepixel_gpu_wsl.ps1 `
  -Project elmer_project_gpu_hybrid_prism_phase19.json `
  -Case case_p19_gpu_amgx_phase19_time5us_smoke_7step
```

本番の比較区間（20.020 ms のパルス後 100 us、5 us 以下の時間刻み）は次のケースである。

```powershell
.scriptsun_singlepixel_gpu_wsl.ps1 `
  -Project elmer_project_gpu_hybrid_prism_phase19.json `
  -Case case_p19_gpu_amgx_phase19_time5us `
  -AmgxConfig config\amgx\tes_fgmres_aggregation_l1_1e-9.json
```

出力は `results/<case>/solver.log` と同ディレクトリの系列CSVに保存される。COMSOL比較は次で作成できる。

```powershell
python scripts/analysis/compare_singlepixel_amgx_comsol.py `
  --elmer results\case_p19_gpu_amgx_phase19_time5us\case_p19_gpu_amgx_phase19_time5us_series.csv `
  --out artifacts\comparison\comsol_gpu_amgx_phase19_time5us_100us `
  --end-us 100 `
  --solver-label "AMGX / RTX 3060 Ti / Phase19 hybrid-prism"
```

## 判定

GPUバイナリの MUMPS（AMGXなし）は同一Phase19メッシュで完走し、CPU/MUMPSと同じ系列を再現した。一方、AMGXは疎行列の収束表示が成功しても、モルタル拘束を含む熱方程式で解の誤差が温度感度 `alpha` に増幅される。production-v2 で高ペナルティ（`1e8`）を使った177ステップ実行は完走したが、COMSOLとの差は最大 **5.314 µA**（ピークの68.3 %）で、比較合格とはしない。

追加診断では、行平衡化とElmer本体の系スケーリングを無効化した `-AmgxConstraintMode no-scaling` は見かけ上完走するものの、初期残差が許容値以下となり反復0回になる。`1e-9` 設定でもパルス後0.05 µs時点で約0.27 µAの電流変化（CPU/MUMPSは約0.0004 µA）を生じ、残差判定が物理誤差を検出できないことが確認できた。未スケール系は比較用途には不適である。

収束判定を `5e-12` に厳しくした通常スケール系も2000反復で残差 `2.56e-10` に留まり、反復回数を増やすだけでは解決しない。

既存の `slave` 縮約モードも10ステップで2000反復・残差 `4.68e-9` に停滞した。既存のSchurハイブリッドは1ステップで外側反復が残差 `6.5e-3` 付近から進まず、現状のまま本番計算へは使わない。

短時間診断のプロットと数値は `artifacts/comparison/comsol_gpu_amgx_prod_v2_noscaling_smoke10_0p05us/` に保存した。

したがって現状は「GPU実行・ログ・比較プロットまで」は完成しているが、「COMSOLと同じ物理波形を保つGPU線形ソルバ」には未到達である。原因はメッシュや時間刻みではなく、AMGXへ渡す前のモルタル拘束（消去／ペナルティ）の数値安定性である。AMGXを比較用途へ昇格するには、拘束付き系のSchur補正またはCPU側拘束解法とのハイブリッド化を追加し、まず同一メッシュのCPU/MUMPSとの差を閾値（基線1%、波形2%）内にする必要がある。
