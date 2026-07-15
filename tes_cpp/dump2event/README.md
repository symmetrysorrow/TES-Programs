# tes_cpp dump2event

`dumpall.dat` をイベント単位の `event.json` に変換する、依存なしの C++17 ライブラリと `dump2event` CLI です。Pulse 計算・Eigen・Visual Studio 固有の設定には依存しません。

## CMake での利用

```powershell
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --config Release
.\build\Release\dump2event.exe dumpall.dat event.json --input-energy 1.0
```

`--input-energy` の単位は dumpall.dat 内のエネルギーと同じ単位です。既存アプリで `InputPara.E / 1000` を渡しているため、通常は MeV です。

すべてのイベントを保存する場合は次のようにします。

```powershell
dump2event dumpall.dat event.json --input-energy 1.0 --save-all --full-energy-list FullEnergyList.dat
```

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
