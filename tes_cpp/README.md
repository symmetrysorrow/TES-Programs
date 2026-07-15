# tes_cpp

PHITS の `dumpall.dat` 変換とパルス生成を同じ C++／Python パッケージで提供します。CLI は `dump2event` と `posi2pulse` です。

```python
from tes_cpp import posi2pulse

# Python のリストとして受け取る
pulses = posi2pulse("input.json", [1, 5, 10])

# JSON ファイルに保存する（戻り値は Path）
path = posi2pulse("input.json", [1, 5, 10], output_path="pulses.json")
```

```python
from tes_cpp import dump2event

dump2event("dumpall.dat", "event.json", input_energy=1.0)
```

position は既存コードと同じ、1 始まりの absorber block 番号です。JSON 出力は、再現用の
`input`、共通の `time`、位置番号をキーとする `pulses` で構成されます。

```json
{
  "input": { "...": "input.json の内容" },
  "time": [0.0, "..."],
  "pulses": {
    "1": { "ch0": ["..."], "ch1": ["..."] }
  }
}
```

## `pulse_model.py` との一致確認

`simu_py/input.json` を両実装に与え、`position` に含まれる全位置の
CH0/CH1 と時刻配列を比較します。先に `posi2pulse` をビルドしてから実行します。

```powershell
python tes_cpp/tests/compare_with_pulse_model.py --executable build/tes_cpp/posi2pulse/Release/posi2pulse.exe
```

既定の許容差は相対誤差 `1e-8`、絶対誤差 `1e-16` です。別の同形式の
入力は `--input path/to/input.json`、別のビルド成果物は `--executable` で指定できます。

```powershell
# TES-Programs のルートから実行する。build/ は配布物に含めない。
cmake -S tes_cpp -B build/tes_cpp -DCMAKE_BUILD_TYPE=Release
cmake --build build/tes_cpp --config Release

# Python パッケージとして使う場合
python -m pip install ./tes_cpp
```

ビルド時には CMake が Eigen 3.4、nlohmann/json 3.11 を自動取得します。ネットワークを使わない環境では、事前に CMake の `find_package` で検出できる場所に両依存を導入してください。

`dump2json/` と `posi2pulse/` は兄弟の CMake コンポーネントです。ルートの
`tes_cpp/CMakeLists.txt` は両方をまとめてビルドし、各サブディレクトリは単独でも
構成できます。

生成された `build/`、`tes_cpp/build-*` と `CMakeCache.txt` はマシン固有の
絶対パスを含むため、配布物には含めず、展開先で上記の構成コマンドを実行します。
