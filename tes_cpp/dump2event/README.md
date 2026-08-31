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

全エネルギー判定は、履歴中にリーク終了がないことを主条件にします。`E_deposit` は各履歴で、入射エネルギーから `NCOL=12` のリークエネルギーを差し引いた値に閉じるよう補正されます。これにより、負の局所沈着値を出力せず、リークのない履歴では沈着エネルギーの総和が入射エネルギーになります。

すべてのイベントを保存する場合は次のようにします。

```powershell
dump2event dumpall.dat event.json --input-energy 1.0 --save-all --full-energy-list FullEnergyList.dat
```

初回解析で全吸収イベントだけを保存し、リーク (`NCOL=12`) を検出した履歴をその場で破棄する場合は次のようにします。

```powershell
dump2event dumpall.dat event.h5 --input-energy 0.663 --full-energy-only
```

このモードでは、リークが確定した履歴について残りの座標・エネルギー配列を作りません。ただし、全吸収かどうかを判定するため、`dumpall.dat` の入力行を履歴の終端まで読む処理は残ります。

履歴ごとの粒子終端・反応内訳を確認する場合は、`--history-summary history_summary.csv` を追加します。これは、リークした光子が一次光子か二次光子か、履歴に反応や二次粒子があったか、沈着エネルギーがいくらかをCSVに出力します。

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
