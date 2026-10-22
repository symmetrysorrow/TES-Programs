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
  --physics-config pb_config.json
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

## Pb / Bi material preparation

The converter now recognizes **Pb and Bi** as the only selectable material
identities. Pb remains the only material with approved PHITS same-run internal
score verification. Bi has an explicit research-only material/medium mapping,
source guard and HATCH SHA-256 input validation, but **normal Bi HDF5 v2
publication is blocked** until independent Bi `[t-deposit]` oracle matching.
Do not substitute Pb EDedx/PDedx or density for Bi; do not run expensive new
10,000-history simulations to bypass this gate. See `BI_VALIDATION_PLAN.md`.

### Pb HATCH table included by default (after source-data import)

The preferred setup is to keep the **exact PHITS 3.341/Pb HATCH CSV** in
`tes_cpp/dump2event/data/pb_egs5_phits3341.csv`, install it with the package,
and automatically verify its SHA-256 before using it. The validated historical
file is exactly 130,754 bytes and has the SHA-256:

```text
20fd86eb0b57483bd19f20a8388d2c1ecb0431a02fa99f0933e2c3c4748514e9
```

`--stopping-table <path>` remains supported for private local tables and
legacy workflows, but when omitted the native CLI delegates to the Python
driver to resolve the verified bundled resource. No guessing or network fetch
occurs: missing/changed bundled data cause a hard error. The table bytes have
**not** yet been committed because this repository is public, the actual
historical CSV is available only on a private local Windows installation,
and permission to publicly redistribute PHITS-generated coefficients has
not yet been checked. See `data/README.md` and
`tools/import_verified_pb_table.py`. Until the verified data file is imported,
use the previous explicit `--stopping-table` argument.

#### Private Pb table setup

The exact historical HATCH CSV remains private because its redistribution
rights have not been established. Set `TES_CPP_PB_HATCH_TABLE` once to the
verified local file; the normal CLI then needs no table path:

```powershell
$env:TES_CPP_PB_HATCH_TABLE = 'D:\private\egs5_stopping_power_coefficients.csv'
dump2event input.dumpall output.h5 --physics-config pb_config.json
```

The Python API and PoST UI use the same resolver. Explicit table paths are
also accepted only when their bytes match the pinned SHA-256 above. Missing
or modified coefficients fail before conversion; HDF5 v2 provenance records
the SHA-256 of the actual table used.
