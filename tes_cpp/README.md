# tes_cpp

PHITS 3.341のPb／EGS5物理deposit再構成と、独立したTESパルス関連ツールを
提供するC++／Pythonパッケージです。

## dump2event（HDF5 v2専用）

```python
from tes_cpp import dump2event

result = dump2event(
    "dumpall.dat", "event_pb_egs5_v2.h5",
    physics_config="pb_config.json",
    stopping_table="pb_hatch_private.csv",
)
print(result.output_path)
```

古い反応差分式／総和補正を使う変換器は削除しました。
現在はPb、EGS5、単位重み等の検証済み条件のみサポートし、
条件に合わなければエラーです。`Pulse_Ms`はHDF5 v2に対応しますが、SEGMENTのエネルギー付与位置は線分内一様dE/dsという近似モデルで分配します。
検証済み条件・入力設定は`dump2event/README.md`、
HDF5構造は`dump2event/HDF5_V2_SCHEMA.md`を参照してください。

## posi2pulse

```python
from tes_cpp import posi2pulse
pulses = posi2pulse("input.json", [1, 5, 10])
# 現在のパルスAPI。新しいHDF5 v2の熱配分は未対応です。
```

### ビルド

```powershell
cmake -S tes_cpp -B tes_cpp/build -DCMAKE_BUILD_TYPE=Release
cmake --build tes_cpp/build --config Release
python -m pip install .
```

外部のPHITS配布ソース・private停止能係数表や生Dumpallはパッケージへ
組み込まず、実行時にユーザーが指定します。
