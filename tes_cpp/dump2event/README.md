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

tes_cpp/build/dump2event/Release/dump2event.exe input.dumpall output_v2.h5
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
**`Pulse_Ms()`はHDF5 v2のPOINT/SEGMENTを読み込めます。**
POINTはPHITSスコア座標で、SEGMENTは始点〜終点の線分を各absorber blockと
交差させ、交差した線分長に比例する割合でエネルギーを分配します。
これは線分内の停止能(dE/ds)が一様という**明示的な近似**であり、
EGS5の真の空間熱分布を検証したことにはなりません。
不完全/未解決イベント、Pb以外の物質、座標範囲外depositは拒否します。
全吸収判定がundeterminedでも、既知の付与エネルギーパルスは計算できます。

各位置のイベントファイルは`<event-root>/<posi>/event.h5`または
`<event-root>/<posi>/event_pb_egs5_v2.h5`に置いてください。1位置のみなら
`--event-root`にHDF5ファイル自体を指定可能です。複数位置を一度に
dump2event→Pulse_Msへ渡す場合は、先に位置ごとに変換してください。
出力`pulse_MS.h5`の属性に空間配分モデルを記録します。

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

### Bundled Pb HATCH table

The preferred setup is to keep the **exact PHITS 3.341/Pb HATCH CSV** in
`tes_cpp/dump2event/data/pb_egs5_phits3341.csv`; package builds install it and
automatically verify its SHA-256 before using it. The validated historical
file is exactly 130,754 bytes and has the SHA-256:

```text
20fd86eb0b57483bd19f20a8388d2c1ecb0431a02fa99f0933e2c3c4748514e9
```

`--stopping-table <path>` remains available for diagnostics, but is not needed
for normal use. The native CLI delegates to the Python driver, which resolves
the verified bundled resource. No guessing or network fetch occurs: missing
or changed bundled data cause a hard error. Permission from a PHITS
administrator was obtained for individual-use scope; the GitHub repository is
private and currently has only the owner as a collaborator. Do not publish the
CSV, mirror it, or expand access beyond that scope. See `data/README.md` and
`tools/import_verified_pb_table.py`.

#### Private Pb table setup

For installations that do not include the bundled resource, set
`TES_CPP_PB_HATCH_TABLE` to the verified local file. On this branch, the
bundled resource takes precedence and the environment variable is unnecessary:

```powershell
$env:TES_CPP_PB_HATCH_TABLE = 'D:\private\egs5_stopping_power_coefficients.csv'
dump2event input.dumpall output.h5
```

The Python API and PoST UI use the same resolver. Explicit table paths are
also accepted only when their bytes match the pinned SHA-256 above. Missing
or modified coefficients fail before conversion; HDF5 v2 provenance records
the SHA-256 of the actual table used.


## Automatic Pb run association

The normal two-path command and Python call need no physics JSON or table path:

```powershell
dump2event input.dumpall output.h5
```

Automatic mode accepts exactly one sibling `*.inp` whose `[Parameters] file(15)`
resolves to that Dumpall and whose configured `[Parameters] file(6)` log is
present and reports `Version = 3.341`. The input must pass the existing Pb
source, geometry, material, cutoff, tally, seed, and unit-weight checks.
Relative run output paths follow the co-located convention: PHITS input, log,
and Dumpall share a directory. Missing or ambiguous pairs fail with an error.
The generated effective contract, its SHA-256, input/log hashes, source seed,
and bundled table hash are stored in HDF5 v2 provenance.

`--physics-config` and `--stopping-table` remain available for diagnostics.
Explicit contracts still must identify the uniquely associated input and pass
the same source and physics validators. PoST asks for Dumpall and HDF5 output;
it only asks for a private CSV if the verified bundled table is unavailable.
