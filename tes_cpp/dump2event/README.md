# dump2event — Pb EGS5 physical deposition → HDF5 v2

**現在の唯一の変換方式は、PHITS 3.341／Pb／EGS5のスコア再構成です。**
旧式の親エネルギー差・負値補正・全履歴エネルギー正規化方式、JSON出力、
HDF5 v1の新規変換出力、`--input-energy`、`--save-all`、
`--full-energy-only`、`--normalize-deposits`、旧reaction-auditオプションは
変換器から削除しました。旧引数の指定はエラーになります。

## 対応する物理条件

PHITS 3.341、Pb region 1000、Dumpall material ID 2、EGS5 medium 2、
1 keV電子・陽電子cutoff、単位重み、同一runに対応する
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

Automatic mode accepts exactly one sibling `*.inp` whose explicit or
PHITS-reported effective `[Parameters] file(15)` resolves to that Dumpall and
whose configured `file(6)` log (default `phits.out`) reports `Version = 3.341`.
PHITS defaults are resolved before validation. For EGS5, omitted
`emin(12:14)` values use the PHITS defaults for `negs=1/2`. PHITS 3.341
also aligns `emin(12)` to `emin(13)` when they differ; the resolver records
that effective value and its source. In the `dump_662` input, explicit
`emin(13)=0.001` causes omitted `emin(12)` to be adjusted to 0.001 MeV, which
matches the validated electron cutoff. If both charged cutoffs are omitted,
their 0.1 MeV defaults remain outside the current validated scope and are
rejected. See the [PHITS EGS5 parameter defaults](https://phits.jaea.go.jp/manual/PHITS-en/chapters/sections/parameters/parameters_13.html).
The input must pass Pb source type/energy, geometry, material, charged-particle transport mode,
cutoff, and unit-weight checks; `[t-deposit]`, an explicit `bitrseed`,
`maxcas`, and `maxbch=1` are no longer mandatory.
Relative run output paths follow the co-located convention: PHITS input, log,
and Dumpall share a directory. Missing or ambiguous pairs fail with an error.
The generated effective contract, its SHA-256, input/log hashes, source seed,
and bundled table hash are stored in HDF5 v2 provenance.

`--physics-config` and `--stopping-table` remain available for diagnostics.
Explicit contracts still must identify the uniquely associated input and pass
the same source and physics validators. PoST asks for Dumpall and HDF5 output;
it only asks for a private CSV if the verified bundled table is unavailable.


## Gamma-ray source energy scope: (0, 2] MeV

The source-input validator accepts **any finite, positive incident photon energy
up to and including 2.0 MeV**. A previously hardcoded allowlist
(Co-60 1.173/1.332 MeV, 0.511 MeV, 2.0 MeV) no longer blocks intermediate
energies such as 0.662 or 0.663 MeV. These are **input-syntax** acceptances,
not proof of new PHITS same-run deposit accuracy at each energy.

For generic discrete photon lines, require `[Source] proj=photon`,
`e-type=8`, integer `ne` equal to the actual number of energy lines,
`totfact=1`, `wgt=1` if provided, and intensity **1 on every line**.
Non-finite, zero, negative, and >2-MeV energies fail closed. The exact
historical two-line Co-60 syntax is retained for backward compatibility.
Source energies and modes remain part of the resolved physics provenance.

`[t-deposit]` is a separate **deposit-energy histogram** setting, not the
input gamma energy and not the EGS5 electron/positron stopping-power table.
The validated legacy `emax=1.5` histogram stays acceptable, including for
2-MeV input, to avoid invalidating past runs; it does **not** cover deposits
above 1.5 MeV in the histogram. When choosing a new histogram bin range,
set `emax` at or above the largest source photon energy (e.g.
`emax=2.0` for 2-MeV input). The `[t-deposit]` region, output, unit,
and other established strict gates are unchanged.

The existing source-input pairing (same-folder effective `file(15)` Dumpall and
`file(6)` PHITS 3.341 log), explicit 64-bit `bitrseed`, Pb material
and density, 1-keV charged cutoffs, SHA-pinned Pb HATCH table,
unit-weight requirements, and unsupported-process `UNRESOLVED` handling
continue to apply. In particular a `dump_662` run with no matching
`[t-deposit]` or seed **still fails**, regardless of its accepted energy.

Physical scoring equivalence was demonstrated for previously validated
inputs, not automatically for all gamma energies in this numeric interval.


## Minimal required PHITS input for V2 conversion

The Dumpall-derived energy calculation does **not** consume a PHITS
`[t-deposit]` file, any `[t-deposit]` histogram, the PHITS random
`bitrseed`, or the requested number of histories (`maxcas`). Those are
no longer mandatory declarations. `maxbch` is metadata only; the event
exporter preserves Dumpall event order; NOCAS identifiers must be positive and
unique, but they do not need to be numerically increasing. Duplicate
identifiers are rejected rather than guessed to represent different batches.

The `[t-deposit]` section may be absent or contain a different histogram
specification: its presence and binning do not alter the reconstructed
physical ledger. When present, the provenance reports
`t_deposit_status=configured_not_used`; when absent, `absent`.
No PHITS tally results are silently used or inferred.

An absent `bitrseed` is recorded as
`source_seed_status=not_explicit_in_input` and
`source_seed_bits=null` in the generated contract; no seed is fabricated.
When a seed is explicitly supplied, its 64-bit binary format is checked.
The source-input and Dumpall byte SHA-256 hashes remain mandatory provenance.
If `maxcas` and `maxbch` are present, their positive integer values are
recorded; an absent value is `null`.

Optional `[Source]` declarations `totfact` and `ne` may be omitted,
but explicit non-unit `totfact` and contradictory `ne` still fail.
An unambiguous discrete photon source (`proj=photon`, `e-type=8`,
positive finite line energy `<=2 MeV`, unit line intensities), the exact
Pb geometry/material, PHITS 3.341 run version, EGS5 electron transport
mode, 1 keV relevant charged-particle cutoffs, and **unit weights on the
actual Dumpall records** remain required. The PHITS source-input pairing
still requires a unique `file(15)` Dumpall target and a `file(6)` run log
confirming PHITS 3.341.

Removing a tally or source metadata requirement does **not** extend the
physics process validation: unknown additions remain UNRESOLVED and prevent
publishing a falsely complete HDF5 v2. The existing same-run PHITS scoring
comparison only supports the previously checked runs, not arbitrary
<=2-MeV photon source geometries or physics processes.
