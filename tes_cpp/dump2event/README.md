# dump2event — Pb EGS5 physical deposition → HDF5 v2

**現在の唯一の変換方式は、PHITS 3.341／Pb／EGS5のスコア再構成です。**
旧式の親エネルギー差・負値補正・全履歴エネルギー正規化方式、JSON出力、
HDF5 v1の新規変換出力、`--input-energy`、`--save-all`、
`--full-energy-only`、`--normalize-deposits`、旧reaction-auditオプションは
変換器から削除しました。旧引数の指定はエラーになります。

## 対応する物理条件

PHITS 3.341、Pb region 1000、Dumpall material ID 2、EGS5 medium 2、
1 keV電子・陽電子cutoff、単位重み、対象`[t-deposit]`設定、同一runに対応する
明示的Pb HATCH停止能テーブル。これ以外は検証済みとみなさず、fail-closedで
拒否します。検証できた2 MeV・1,000履歴では同一runの内部加算
**5,978,912/5,978,912**件が一致しました（最大差約4.68e-14 MeV）。
これは任意のPHITS入力・全吸収・SEGMENT内の熱分布の保証ではありません。

## ビルド・実行

```powershell
cmake -S tes_cpp -B tes_cpp/build -DCMAKE_BUILD_TYPE=Release
cmake --build tes_cpp/build --config Release

tes_cpp/build/dump2event/Release/dump2event.exe input.dumpall output_v2.h5 `
  --physics-config pb_config.json --stopping-table pb_hatch_private.csv
```

環境によって実行ファイルの位置が異なります。
従来の明示オプション `--physics pb-egs5 --hdf5-version 2` も移行用に
受け付けますが、変換方式の選択肢は一つだけです。

設定JSONは `validation/pb_physics_cli_integration/pb_config.example.json` を
参考にしてください。JSONの`phits_input`は実在するPHITS入力を指し、
停止能CSVはユーザーが保有するprivateな正しい表でなければなりません。
PHITSソース、元Dumpall、HATCH実測表、内部traceは配布しません。

通常Python API：

```python
from tes_cpp import dump2event

result = dump2event(
    "dumpall.dat", "event_pb_egs5_v2.h5",
    physics_config="pb_config.json",
    stopping_table="pb_hatch_private.csv",
)
print(result.output_path)
```

従来の`input_energy`等は拒否します。出力先が既存の場合も拒否します。
CLIはnative ledger・transport exporter・Pythonドライバーを使用するため、
Pythonパッケージから利用する場合はリポジトリルートの`python -m pip install .`
で必要な補助実行ファイルをインストールしてください。

## UI・下流の制約

`PoST_Simulations/PoST_Simulation.py --interactive` の
`Dump2Event`を選ぶと、物理条件・停止能表・入出力を明示してHDF5 v2を生成します。
旧変換メニューはなくなりました。

HDF5 v2はPOINT／SEGMENT／UNRESOLVED、PHITS代表スコア点・区間端点、
履歴・粒子ID・provenanceを区別して保存します。
**`Pulse_Ms()`は未変更の旧点付与モデルであり、HDF5 v2を明示拒否します。**
新しい`event2pulse`がPOINT/SEGMENTの空間熱配分を検証して実装するまで
自動接続しません。全吸収判定も証明されていない場合はundeterminedです。

スキーマは`HDF5_V2_SCHEMA.md`、物理検証・適用範囲は
`validation/dead_secondary_threshold_and_gamma_10k/FINAL_DECISION.md`を参照してください。
