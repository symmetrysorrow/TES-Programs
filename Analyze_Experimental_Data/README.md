# TES Programs / Analysis

## 起動

プロジェクトのルートで次を実行します。

```powershell
python main.py
```

`main.py` がデータフォルダの種類を判定し、Pulse / IV / RT の適切な処理へ振り分けます。

## フォルダ構成

- `main.py`: アプリケーションの唯一の起動入口
- `tes_analysis/`: アプリ本体
  - `cli.py`: 対話メニューと実行フロー
  - `operations.py`: Pulse、Noise、校正、RTの処理
  - `analysis_utils.py`: 共通の読み込み、信号処理、解析関数
  - `dispatch.py`: データ種別の判定と振り分け
  - `prompts.py`: 対話プロンプト
  - `iv.py`, `rt.py`: IV / RT固有処理
- `scripts/`: 個別に実行する補助スクリプト
- `tests/`: 比較用スクリプトと測定データ

補助スクリプトは必要な場合だけ、ルートから次のように実行します。

```powershell
python -m scripts.list_pulses
python -m scripts.inspect_noise
```

## PulseConfig.json の解析パラメータ

`Analysis` ブロックで使えるキーのうち、エネルギー分解能に効くもの。

| キー | 既定 | 意味 |
| --- | --- | --- |
| `CutoffFrequency` | 必須 | Bessel ローパスのカットオフ [Hz] |
| `BaseStart` | 未設定 | ベースライン区間の開始を `PreSample - BaseStart` で指定（Getpara の `base_x`） |
| `BaseWidth` | 未設定 | ベースライン区間の長さ（Getpara の `base_w`） |
| `PeakSearchSample` | 未設定 | 波高探索を `[PreSample, PreSample+PeakSearchSample)` に限定 |
| `PeakAveragePreSample` / `PeakAveragePostSample` | 必須 | 波高平均の窓 |

`BaseStart` / `BaseWidth` を省略すると従来どおり `pulse[0:PreSample]` 全体の
平均をベースラインにします。Getpara と同じ区間にそろえたい場合は、
`presamples=1000` に対して次のように書きます。

```json
"BaseStart": 1000,
"BaseWidth": 500
```

ベースラインはゲイン補正 (`TempCalib`) の説明変数そのものなので、区間の
取り方が `Base` のばらつき、ひいては補正後の分解能に直接効きます。

## 最適フィルタ

方式は理論的な最適フィルタ `H(f) ∝ S*(f)/PSD(f)` の一本だけです。実データで
`1/ASD` 重み付けやテンプレートへの Bessel を含む方式と比較した結果、こちらが
明確に最良だったため、他の方式は削除しました。

- 重み: `1/PSD`（`PSD = ASD**2`）
- テンプレートには Bessel を掛けない（各パルス側に掛かっているので、
  掛けると帯域制限が二重になる）
- `Σ|S|²/PSD` で規格化。平均パルスと同じ波形に対して推定量がちょうど 1 に
  なり、そこへ平均パルスの波高を掛けて `Peak` と同じスケールに戻す

`modelnoise.txt` に保存されるのは **ASD**（振幅密度, pA/√Hz）です。規格化の
おかげで、ノイズモデル全体に掛かる定数倍（eta や FFT 正規化の流儀、ASD と PSD
の取り違え）は推定量に影響しません。効くのは PSD の *形* だけです。

出力される波高推定量:

| カラム | 内容 |
| --- | --- |
| `Peak` | Bessel 後の波高平均 |
| `PeakOpt` | 最適フィルタ出力。`Peak` と同じスケール |
| `PeakTemp` / `PeakOptTemp` | 上記に `TempCalib`（ベースライン依存ゲイン補正）を掛けたもの |

`Compare Estimators` メニューで、これらの mean / std / FWHM / FWHM÷mean を
一覧表示し、`Base` vs 各推定量、`Decay` vs `PeakOpt`、ヒストグラムをまとめて
確認できます。`operations.CompareChain()` は Average Pulse → Noise spectrum →
Template → 各イベント出力 → 補正前 FWHM → 補正後 FWHM の順に、どの段階で
差が付いているかを表示します。
