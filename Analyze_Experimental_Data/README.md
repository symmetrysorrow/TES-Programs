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
