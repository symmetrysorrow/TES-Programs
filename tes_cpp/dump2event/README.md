# tes_cpp dump2event

`dumpall.dat` をイベント単位の `event.json` に変換する、依存なしの C++17 ライブラリと `dump2event` CLI です。Pulse 計算・Eigen・Visual Studio 固有の設定には依存しません。

## CMake での利用

```powershell
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --config Release
.\build\Release\dump2event.exe dumpall.dat event.json --input-energy 1.0
```

`--input-energy` の単位は dumpall.dat 内のエネルギーと同じ単位です。既存アプリで `InputPara.E / 1000` を渡しているため、通常は MeV です。

通常運用で指定するのは `--input-energy`、`--save-all` または `--full-energy-only`、必要なら `--full-energy-list` です。`--history-summary`、`--max-histories`、`--summary-only` は診断・開発用です。`PoST_Simulation.py` から実行する場合、これらのCLIフラグを直接指定する必要はありません。

全エネルギー判定は、履歴中にリーク終了がないことを主条件にします。**デフォルトではエネルギー総和を入射エネルギーに強制一致させません。** 反応点の `E_deposit` は、親トラックの直前点エネルギー `E` から輸送される二次粒子の運動エネルギー合計を引いた値として記録します（PHITSの `QCLUSTS(7)`、`JCLUSTS(4)` と輸送閾値で集計）。二次粒子が0個の場合は親トラックの `E` を局所付与量に記録します。負の付与量は隠さず残し、元データの矛盾として確認してください。

従来の負値切り捨てと `input_energy - leaked_energy` への全履歴比例正規化は、再現用の **`--normalize-deposits`** 指定時のみ実行します。このモードでも `--history-summary` の `raw_total_deposit`、`negative_deposits`、`energy_consistent` は補正前の値で評価します。 `energy_consistent` は `raw_total_deposit` と `source_energy - leaked_energy` の差に基づきます（一次エネルギーが読めない場合は `input_energy` を使用）。リークを検出して早期破棄した履歴の診断値は未完了です（`deposit_records_complete=0`）。エネルギー収支や二次粒子の解釈は、実際のPHITSダンプと[t-deposit]による独立照合が必要です。

すべてのイベントを保存する場合は次のようにします。

```powershell
dump2event dumpall.dat event.json --input-energy 1.0 --save-all --full-energy-list FullEnergyList.dat
```

初回解析で全吸収イベントだけを保存し、リーク (`NCOL=12`) を検出した履歴をその場で破棄する場合は次のようにします。

```powershell
dump2event dumpall.dat event.h5 --input-energy 0.663 --full-energy-only
```

このモードでは、リークが確定した履歴について残りの座標・エネルギー配列を作りません。ただし、全吸収かどうかを判定するため、`dumpall.dat` の入力行を履歴の終端まで読む処理は残ります。

履歴ごとの粒子終端・反応内訳を確認する場合は、`--history-summary history_summary.csv` を追加します。これは、リークした光子が一次光子か二次光子か、履歴に反応や二次粒子があったか、沈着エネルギーがいくらかをCSVに出力します。追加列 `raw_total_deposit`、`negative_deposits`、`deposit_records_complete`、`normalization_factor` を用いると、正規化前のエネルギー収支とデータ欠落を検証できます。さらに `source_energy` と `source_energy_available` は `NCOL=4` の一次粒子から読み取った履歴ごとの入射エネルギーです。異なるγ線が混在する場合、`energy_consistent` は可能な限りこの実測値で計算します（`source_energy_available=0` ならCLIの `--input-energy` にフォールバック）。 `energy_balance_reference` は利用した値です。現状は各履歴に単一の一次源粒子があるケースを想定しており、同じ履歴で複数の独立した源粒子が放出される設定は別途検証が必要です。

先頭の履歴だけを調べる場合は、`--summary-only --max-histories 1000` を追加するとイベントファイルを書かずに指定件数で停止できます。

## Python での利用

wheel を作成・導入するとネイティブ CLI も同梱されます。

```powershell
python -m pip install .
```

```python
from tes_cpp_dump2event import convert

result = convert("dumpall.dat", "event.json", input_energy=1.0, save_all=True)
print(result.full_energy_event_ids)
```

開発中に CMake でビルドした実行ファイルを使う場合は、`TES_CPP_DUMP2EVENT_EXECUTABLE` にそのパスを設定してください。
