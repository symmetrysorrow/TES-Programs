# TES ID と `show` CLI 再設計 実装計画

日付: 2026-09-29  
状態: Planned

## 目的

現在の schema-v3 project/model TOML では、TES のユーザー向け ID が model の内部構造から自動生成されている。

例:

- `pixel.tes`
- `west.a`
- `west.b`
- `east.a`
- `east.b`

これらは `layout.groups[].id` と component 内のローカル TES ID を結合した値であり、model の組み立て方を知らないユーザーには意味が分かりにくい。また、`python main.py show <project>` は現在、project/model を解決した結果の一部しか表示せず、TOML を直接見ることとの差が小さい。

本変更では次の2点をまとめて改善する。

1. **TES の公開 ID と model 内部 ID を分離する**
2. **`show` を「resolved configuration viewer」として再設計する**

目標は、project を編集するユーザーが model 内部の group/component 構造を意識せずに TES を参照でき、`show` を見れば「この計算が最終的にどう解釈されるか」を把握できる状態にすることである。

---

## 設計方針

### 1. 公開 TES ID と内部 ID を分離する

component 内の TES ID と layout group ID は、model を構成するための内部識別子として残す。

例:

```toml
[[components.tes_pair.tes]]
id = "a"
position = ["-2[mm]", "0[mm]"]

[[components.tes_pair.tes]]
id = "b"
position = ["2[mm]", "0[mm]"]
```

一方、project から参照する TES ID は layout 側で明示的に割り当てる。

```toml
[[layout.groups]]
id = "west"
component = "tes_pair"
origin = ["-5[mm]", "0[mm]"]
tes_ids = { a = "tes_1", b = "tes_2" }

[[layout.groups]]
id = "east"
component = "tes_pair"
origin = ["5[mm]", "0[mm]"]
tes_ids = { a = "tes_3", b = "tes_4" }
```

これにより、

- `west`: layout 上の group ID
- `a`: component 内のローカル TES ID
- `tes_1`: project から参照する公開 TES ID

という役割分担になる。

公開 ID は位置や並び順から自動生成しない。geometry の変更で TES の順序が変わっても identity が変化しないよう、model 作者が明示的に定義する。

### 2. `tes_ids` は mapping とする

次のような配列形式は採用しない。

```toml
tes_ids = ["tes_1", "tes_2"]
```

component 内の TES 順序に意味が発生し、並べ替えによる誤対応が起きるためである。

必ずローカル ID と公開 ID の対応を明示する。

```toml
tes_ids = { a = "tes_1", b = "tes_2" }
```

### 3. model version を 3 に上げる

TES identity の決定規則が変わるため、model format は `model_version = 3` とする。

互換性は次のように扱う。

| model version | TES ID の決定方法 |
| --- | --- |
| 2 | 従来通り `<group>.<local_id>` |
| 3 | `layout.groups[].tes_ids` の mapping |
 
リポジトリ内の現行 model は v3 へ移行する。v2 サポートは既存入力の読み込み互換性のために残し、新規 model の推奨形式にはしない。

project 側の構造は `circuits.<TES ID>` のままであり、conceptual schema は変わらないため、project の `schema_version = 3` は維持する。

---

## `scenario_config.py` の変更

### TES 展開

v3 model では、group ごとに component TES を展開するとき、公開 ID を `tes_ids` から取得する。

概念的には次の処理になる。

```text
component local id "a"
    ↓
group.tes_ids["a"]
    ↓
public id "tes_1"
```

### 検証

v3 model では以下を検証する。

- `tes_ids` が table/mapping である
- `tes_ids` の key が component 内 TES ID と完全一致する
- 必要な local ID の mapping が不足していない
- component に存在しない local ID が mapping に含まれていない
- 公開 TES ID が `ID_RE` に適合する
- 公開 TES ID が全 group を通して一意である
- project の `circuits` key が公開 TES ID と完全一致する

エラーは構造を特定できる内容にする。

例:

```text
group west: missing public TES ID for component TES 'b'
duplicate public TES ID: tes_2
unknown circuit TES IDs: ['foo']; available IDs: ['tes_1', 'tes_2', ...]
```

### scenario 内部表現

resolved scenario には公開 ID だけでなく、model 内部の provenance も保持する。

例:

```python
{
    "id": "tes_1",
    "group_id": "west",
    "component": "tes_pair",
    "local_id": "a",
    "x": -0.007,
    "y": 0.0,
    "circuit": {...},
    "sources": {...},
}
```

通常のユーザー向け出力では `id` のみを使い、group/component/local ID は verbose/debug 表示に利用する。

---

## model TOML の移行

### single pixel

現在:

```toml
[[layout.groups]]
id = "pixel"
component = "tes_single"
origin = ["0[mm]", "0[mm]"]
```

変更後:

```toml
model_version = 3

[[layout.groups]]
id = "pixel"
component = "tes_single"
origin = ["0[mm]", "0[mm]"]
tes_ids = { tes = "tes_1" }
```

### four TES

現在の公開 ID:

```text
west.a
west.b
east.a
east.b
```

変更後:

```text
tes_1
tes_2
tes_3
tes_4
```

対応は model に明示する。

```toml
[[layout.groups]]
id = "west"
component = "tes_pair"
origin = ["-5[mm]", "0[mm]"]
tes_ids = { a = "tes_1", b = "tes_2" }

[[layout.groups]]
id = "east"
component = "tes_pair"
origin = ["5[mm]", "0[mm]"]
tes_ids = { a = "tes_3", b = "tes_4" }
```

番号は位置から自動決定しない。

---

## project TOML の移行

次のような既存記述を、

```toml
[circuits."pixel.tes"]
```

以下へ変更する。

```toml
[circuits.tes_1]
```

4 TES project も同様に、

```toml
[circuits.tes_1]
[circuits.tes_2]
[circuits.tes_3]
[circuits.tes_4]
```

へ移行する。

少なくとも以下を更新する。

- `projects/single_pixel.toml`
- `projects/single_pixel_alpha240.toml`
- `projects/post_four_tes_nominal.toml`
- `projects/post_four_tes_low_bias.toml`

---

## `show` CLI の再設計

### 目的

`show` は TOML の短縮表示ではなく、project と model を読み込んだ後の **resolved configuration** を表示する。

ユーザーが答えを得たい質問は、

> この project を実行したら、どの model / geometry / mesh / TES / 回路値 / pulse / time condition が実際に使われるのか？

である。

### デフォルト表示

表示は概ね以下のセクションに分ける。

```text
Project: post_four_tes_nominal
File:    projects/post_four_tes_nominal.toml
Model:   projects/models/post_four_tes.toml

Geometry
  Absorber: 16 x 1 mm
  TES: 4

Mesh
  Global:   50 um
  TES:       8 um
  Absorber: 35 um
  Layers: stycast=32, Si=4, SiO2=2, SiNx=2

Run
  Pulse:      1332 keV
  Start:      20.02 ms
  Duration:   1 ns
  Sigma:      50 um
  Center:     auto
  Window:     default
  Save field: no

TES
  tes_1  (-7.000, 0.000) mm
    I_bias = 715 uA
    R_sh   = 3.9 mohm
    L_tes  = 12.3 nH
    R_0    = 15.527 mohm
    alpha  = 256.46
    beta   = 5.03
    I_0    = ...

  ...
```

### 表示する値

デフォルト `show` では次を表示対象とする。

- project file
- model file
- absorber の resolved size
- TES 数
- mesh の主要設定
- pulse 条件
- window / timesteps の要約
- save_field
- 各 TES の公開 ID
- resolved position
- `I_bias`
- `R_sh`
- `L_tes`
- `R_0`
- `alpha`
- `beta`
- 導出された `I_0`

solver override が存在する場合は、その project が明示的に変更した主要 solver 条件も表示する。

### 表示しないもの

デフォルトの `show` では model の全実装詳細を列挙しない。

例:

- 全 material property
- 全 geometry parameter
- 内部 template JSON の値
- 生成 case 名の全情報
- provenance/debug 用 source 表示

これらまで出すと model TOML の劣化コピーになるため、詳細は model TOML または verbose/debug 表示に任せる。

---

## `show --verbose`

必要に応じて `--verbose` を追加する。

verbose では、公開 TES ID と model 内部構造の対応を表示する。

例:

```text
tes_1
  layout group: west
  component: tes_pair
  component TES: a
  position: (-7, 0) mm

  I_bias = 715 uA
    source: projects/post_four_tes_nominal.toml
```

既存の `explain()` が持っている source 表示ロジックは、この用途へ統合する。

通常表示では provenance を出さない。

---

## `show` と `run --dry-run` の責務分離

各コマンドの役割を明確にする。

| command | 答える質問 |
| --- | --- |
| `list` | どんな project があるか |
| `show` | この project は最終的にどんな設定になるか |
| `check` | この project は妥当か |
| `mesh` | mesh を準備する |
| `run --dry-run` | 実行すると何をするか |
| `run` | 実行する |

`show` は configuration を表示する。

`run --dry-run` は execution plan を表示する。

dry-run では例えば以下を重視する。

- mesh を再利用するか生成するか
- steady を再利用できるか
- steady -> pulse の実行順
- MPI rank 数
- 実際の生成/計算は行わないこと

mesh hash や内部 case 名などは、通常の `show` の主情報にはしない。

---

## `scenario_cli.py` の整理

現在の `summary()` と `explain()` は責務が重なっているため整理する。

候補:

```python
show_project(scenario, *, verbose=False)
```

または、

```python
show_project(scenario)
show_project_verbose(scenario)
```

`main.py` の mapping も、

```python
{"show": "summary", ...}
```

から、ユーザー向けの明示的な `show` action へ変更する。

内部 CLI に残っている `summary` / `explain` action は、必要性を確認したうえで削除または内部互換 alias とする。

---

## テスト計画

### `tests/test_scenario_config.py`

現在の、

```python
["west.a", "west.b", "east.a", "east.b"]
```

という期待値を、

```python
["tes_1", "tes_2", "tes_3", "tes_4"]
```

へ変更する。

加えて以下をテストする。

1. model v3 の public TES ID が正しく展開される
2. group/local ID が public ID に混ざらない
3. `tes_ids` の不足を reject する
4. `tes_ids` の余分な local ID を reject する
5. public TES ID の重複を reject する
6. 無効な public TES ID を reject する
7. project の未知 TES ID を reject し、available IDs に public ID を表示する
8. circuit 値変更では mesh hash が変わらない
9. public TES ID と circuit key が不一致なら実行前にエラーになる
10. model v2 は従来の `group.local_id` で読み込める
11. model v3 では explicit `tes_ids` を要求する
12. public ID を変えても geometry/mesh の物理内容が同一なら mesh identity が意図せず変化しないことを確認する

最後の項目では、mesh hash の identity に public TES ID を含めるかどうかを明示的に判断する。TES ID は物理 geometry そのものではないため、原則として表示名変更だけで mesh が再生成されない設計を優先する。

### `tests/test_main_cli.py`

現在の、

```python
assert "pixel.tes" in output
```

を、

```python
assert "tes_1" in output
```

へ更新する。

さらに `show` が最低限以下を含むことを確認する。

- model path
- geometry summary
- mesh summary
- public TES ID
- TES position
- `I_bias`
- `R_sh`
- `L_tes`
- `R_0`
- `alpha`
- `beta`
- pulse summary
- window/timestep summary

verbose を実装する場合は、group/component/local ID が verbose でのみ現れることもテストする。

---

## ドキュメント更新

実装後に少なくとも以下を更新する。

- `README.md`
- `projects/README.md`
- `docs/getting_started.md`
- `docs/toml_unification_plan.md`

特に、

> TES ID は model の layout から決まり、show で確認できます

という説明は、次の考え方へ変更する。

> TES の公開 ID は model の layout で明示的に割り当てる。component/group の内部 ID とは独立しており、project は公開 ID のみを参照する。

すべての project 例は `circuits.tes_1` 形式に統一する。

---

## 実装順序

1. `scenario_config.py` に model v3 と `tes_ids` の読み込み・検証を追加する
2. resolved scenario に `id / group_id / component / local_id` を保持する
3. model TOML を `model_version = 3` へ移行する
4. project TOML の circuit key を public TES ID へ移行する
5. `tests/test_scenario_config.py` を更新し、ID validation のテストを追加する
6. `scenario_cli.py` の `show` を resolved configuration viewer として再実装する
7. 必要なら `show --verbose` を追加する
8. `main.py` から新しい `show` action へ接続する
9. `run --dry-run` の表示を execution plan に限定して整理する
10. `tests/test_main_cli.py` を更新する
11. README / getting started / TOML design docs を更新する
12. 全テストを実行する
13. `show -> check -> run --dry-run` の手動 smoke test を行う

---

## 完了条件

以下をすべて満たしたら完了とする。

### ID

- project TOML に `pixel.tes`, `west.a`, `east.b` のような model 内部構造由来 ID が出ない
- project は `tes_1`, `tes_2` など model が明示した公開 ID のみを参照する
- model 内部の group/component/local ID と公開 ID が明確に分離されている
- ID mapping の誤りは実行前に具体的なエラーとして検出される

### show

```console
python main.py show post_four_tes_nominal
```

を見たユーザーが、`west` や `a` の意味を知らなくても、

- どの model を使うか
- TES が何個あるか
- TES がどこにあるか
- 各 TES の最終回路条件
- pulse 条件
- time/window 条件
- mesh の主要条件

を理解できる。

一方で、material 定数や model の全 parameter を確認したい場合は model TOML を直接読む、という役割分担が保たれている。

### CLI の責務

- `show` は resolved configuration
- `check` は validation
- `run --dry-run` は execution plan

という境界がコード・help・docs・テストで一貫している。

---

## 非目標

今回の変更では以下は行わない。

- material/geometry model 全体を `show` に展開する
- 公開 TES ID を座標やソート順から自動生成する
- component/group の仕組み自体を廃止する
- project/model の二層構造を変更する
- circuit override 階層を再導入する
- legacy JSON / legacy TOML の ID 体系を新方式へ移行する

legacy workflow は既存結果の再現経路として現状の名前空間に維持する。
