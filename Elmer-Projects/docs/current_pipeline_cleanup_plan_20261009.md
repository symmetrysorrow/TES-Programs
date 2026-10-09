# 現行 Elmer 計算パイプラインの整理・修正計画

- 作成日: 2026-10-09
- 対象: `Elmer-Projects` の新方式（`projects/*.toml` / `projects/models/*.toml`、schema v3 / model v2）
- 対象モデル: 単ピクセル型、PoST 4TES型
- 調査時の参照コミット: `ec0a418c7b241d9b31dcb983326baeba3aaf8a5c`
- 状態: **ソースコードの静的監査による修正計画**。PoSTの生成済みSIF、実メッシュ、ローカルの改造HeatSolve、実行ログ、COMSOL波形との再比較は未実施。
- この文書のコミットでは**計算コードや設定ファイルは変更しない**。

## 目的と方針

新方式のTOMLは簡素だが、内部で旧JSONテンプレートと多用途SIFビルダーを再利用するため、現在のモデルに不要なレガシー既定値が残っている。特に「完全共形メッシュとMortar」「PoSTと単ピクセルで異なる時間積分」「実行バイナリや再利用キャッシュの不透明さ」は、波形不良の解析を難しくする。

基本方針:

1. **新方式の本番経路を1つの明示的な計算契約にする**。旧非共形・外部回路・診断オプションは `legacy` または別の明示的な実験経路に隔離する。
2. **物理モデルの変更と数値アルゴリズムの変更を分ける**。単ピクセルで検証済みの動作を無条件にPoSTへコピーせず、同一メッシュ・同一状態でA/B比較する。
3. **暗黙フォールバックをなくす**。不明なソルバー名、ビルド不一致、キャッシュ由来不明、非共形接触を検出したら停止する。
4. **本番計算を実行する前に、生成SIF・メッシュ・バイナリ・定常restartを検証する**。

## 修正対象の一覧

| ID | 優先度 | 分類 | 現状／問題 | 推奨アクション |
| --- | --- | --- | --- | --- |
| P0-01 | 最優先 | 数値・界面 | 完全共形を生成してもMortarが既定で有効 | 新方式ではMortar無効を明示、界面節点共有を検証 |
| P0-02 | 最優先 | 数値・時間 | PoSTだけ既定の約20.130001 msで終了、パルス後約110.001 µsしかない。既定BDF1 | 時間窓とBDF次数を独立に明示し、PoSTの長時間応答を計算 |
| P0-03 | 最優先 | 実行環境 | 改造版HeatSolveのロードを新CLIが保証しない | 実行時互換性のpreflightとバージョン指紋を必須化 |
| P0-04 | 高 | 回路連成 | 過去の検証ケースはStep Commitを明示するが新方式は無指定（false） | A/B評価して選定し、採用方式を明示 |
| P1-01 | 高 | 設定検証 | `[solver]` の未知値が拒否されず、誤記したSolver名はUMFPACKへ落ち得る | 新方式の厳密なSolverスキーマを導入 |
| P1-02 | 高 | 再現性 | メッシュと定常restartの再利用条件が緩く、生成器・SIF・実行DLLの変更を全て追跡しない | 内容指紋とマニフェストを検証、履歴不明は再計算 |
| P1-03 | 高 | UDF | UDFコンパイラと実行バイナリの組合せ・同名DLLの検索順が不透明 | 単一のロード元に固定し、実ロードを検証 |
| P1-04 | 中 | メッシュ・パルス | PoSTの中央パルス位置に局所細分化が明示されない | パルス位置の局所細分化と空間収束試験 |
| P1-05 | 中 | モデル寸法 | `Stycast_dx=400 µm` とメッシュ接触直径の既定値 `498 µm` が不一致 | 幾何の唯一の真値を決め、意図を明記 |
| P2-01 | 中 | レガシー継承 | 旧JSONテンプレートを丸ごと内部アダプタにコピー | 新方式専用の最小テンプレートに分離 |
| P2-02 | 中 | 不使用パラメータ | `dabs_dx`、`tes_pitch`、`membrane_dz` などの旧値がモデルとハッシュに残る | 参照先を棚卸しし、不使用値を除外 |
| P2-03 | 中 | ソース生成器 | 旧 `circuit_local` / `circuit_implicit` / `circuit_parallel` / Mortar Block-Schur が現行SIF生成と混在 | レガシー/現行のSIF生成経路を分割 |
| P2-04 | 中 | 診断性 | `save_field=false` が既定で、波形異常時に内部温度・熱流束の解析が困難 | 診断プロファイルで温度場出力と収支を有効化 |

## 詳細と変更案

### P0-01 — 完全共形なのにMortar境界条件が残る

**根拠**

- [`scenario_cli.py:44–79`](../src/support/scenario_cli.py) は新規生成メッシュに `--conformal-tes-stack --conformal-abs` と `ElmerGrid -merge 1e-10` を使用する。
- [`generate_hybrid_prism_geometry.py:737–767`](../src/generate_hybrid_prism_geometry.py) はStycast接触円を吸収体底面へ刻み、同一の面メッシュを写して節点結合する設計。
- [`build_cases.py:1038–1102, 1149–1181, 1695–1701`](../src/support/build_cases.py) では `apply_mortar_bcs=True` が既定。 `conformal_tes_stack=True` が抑制するのはTES–MembraneとStycast–TESのMortarで、Stycast–吸収体は残る。
- 新方式の [`scenario_config.py:294–330`](../src/support/scenario_config.py) は `conformal_tes_stack=True` を旧テンプレートから継承するが、`apply_mortar_bcs=False` は指定しない。

**変更**

- 新方式の定常/過渡に `apply_mortar_bcs=False`、`conformal_tes_stack=True` を明示する（レガシーの既定値は不用意に変更しない）。
- 将来的には `interface_mode="conformal"` を新方式の固定契約とし、Mortarの切替自体を新方式の設定UIから除く。
- ただしSIF上のMortarを外すだけで完了とせず、4TES全接触の節点共有、体積要素の連結性、接触面積、熱流束を実メッシュで照合する。非共形だった場合は**Mortarを外してはいけない**。

**完了条件**: 新方式の単ピクセルとPoSTのSIFに `Mortar BC =` / Galerkin / Plane Projector がなく、`Apply Mortar BCs = False`。全接触の共有節点・連結性が検証される。旧非共形ケースのMortar生成は回帰テストで維持する。

### P0-02 — PoSTだけ時間窓・積分法が変わる

**根拠**

- [`projects/single_pixel.toml`](../projects/single_pixel.toml) と [`projects/single_pixel_alpha240.toml`](../projects/single_pixel_alpha240.toml) は `window = "75ms_bdf2h15"`。
- [`projects/post_four_tes_nominal.toml`](../projects/post_four_tes_nominal.toml) と [`projects/post_four_tes_low_bias.toml`](../projects/post_four_tes_low_bias.toml) は `window`、`timesteps` ともに未指定。
- [`scenario_config.py:313–329`](../src/support/scenario_config.py) の無指定時タイムステップ合計は **20.130001 ms**（パルス開始20.02 msから **110.001 µs** のみ）。`window` を指定した場合だけ `bdf_order=2`、`phase24_bdf2_predictor=True`。
- [`build_cases.py:1483–1497`](../src/support/build_cases.py) では明示しないBDF次数は通常1。
- 従来の[PoST/dual TES記録](dual_tes_plan.md)ではmm単位の伝熱・電流立ち上がりはms域で評価されている。

**変更**

- PoST 2案件に、まず比較実験用として `window = "75ms_bdf2h15"` を設定する。全モデルの波形を観測できる長さを確保する。
- `window`（観測期間/時刻列）と `bdf_order`、BDF2予測子の選択を分離し、無指定によるサイレントな数値方式変更を禁止する。
- 時刻表の検証を追加する（単調、全ステップ正、パルス区間を十分に刻む、予定終了時刻、出力間隔）。
- BDF2の急激なステップ倍率とパルス直後の温度応答についても既存の検証済みスケジュールを参照する。

**完了条件**: 単ピクセルとPoSTで比較可能な時間積分契約を明記し、PoST 4チャンネルのパルス前ベースライン・立ち上がり・ピーク・減衰を同一時間原点で保存する。

### P0-03 — 改造版HeatSolveのロードを保証する

**根拠**

- [`build_cases.py:585–598`](../src/support/build_cases.py) は改造版HeatSolve専用 `TES Inner Circuit Update` をSIFへ書く。
- [`scenario_cli.py:28–40, 169–186`](../src/support/scenario_cli.py) は `tes_parallel_circuit` UDFのコンパイル/配置を行うが、HeatSolveの改造機能の存在検証は見当たらない。
- [`run.py:29–30, 81–132, 837–840`](../src/run.py) のElmerSolver既定値は環境 `ELMER_HOME` に依存し、未設定時は通常のElmer 26.1インストールを指す。SIFヘッダーも `CHECK KEYWORDS Warn`。

**変更**

- 新方式の起動時に「選択されたElmerSolver／実ロードされるHeatSolve／`tes_parallel_circuit`」のパス・SHA256・対応機能を照合する。
- 変更済みHeatSolveが複数回路の `TES Inner Circuit Update` と選定したStep Commit/BDF2処理を提供していなければfail-fast。
- UDFの `elmerf90` は実行するElmerと同じツールチェーンを使い、ロード検証（小さなsmoke runまたは識別子の出力）を追加する。
- `CHECK KEYWORDS Warn` だけに依存してカスタム機能の利用成功を判定しない。

**完了条件**: 標準版Elmerの誤選択は実計算前に検出され、実行マニフェストから採用バイナリを特定できる。

### P0-04 — TES回路のStep Commitを明示する

**根拠**

- [`build_cases.py:510, 597–598, 1623–1648`](../src/support/build_cases.py) の `inner_circuit_step_commit` は未指定時false。
- 過去の[Phase24 Gate 4/5実験](../src/support/run_phase24_gate4_5_nomortar.py) と [Gate 6](../src/support/run_phase24_gate6_gpu_parity.py) は `inner_circuit_step_commit=True` を明示している。
- [`scenario_config.py:310–330`](../src/support/scenario_config.py) は新方式の定常・過渡にこの設定を指定しない。

**変更**

- 同一メッシュ・同一定常状態・同じ時間表でStep Commitのon/offを比較し、TES平均温度・全回路電流・ベースライン連続性を評価する。
- 採用方式を定常/過渡の契約に明示し、falseを偶然の既定値として使わない。
- 未検証のまま機械的にtrueへ変更しない（回路時間離散化は物理波形を変える）。

### P1-01 — `[solver]` を厳密に検証する

**根拠**: [`scenario_config.py:124–128`](../src/support/scenario_config.py) は `solver` が辞書であることしか確認せず、[`build_cases.py:615–754`](../src/support/build_cases.py) の未知 `linear_system` 名は最終 `else` でUMFPACKへフォールバックする。

**変更**: 新方式で使えるSolver種別・各設定の許可キー、数値範囲、バックエンド依存性を宣言して検証。未知キー・未知値は例外。例: `linear_system="mumsp"` はUMFPACKで続行せずエラーとする。

### P1-02 — メッシュ／定常再利用の安全性

**根拠**

- [`scenario_cli.py:44–48, 82–88`](../src/support/scenario_cli.py) は `mesh.names` の存在だけで再利用し、生成時の `PROVENANCE.json` を再照合しない。
- [`run.py:494–519`](../src/run.py) は結果ファイルがあってマニフェストが存在しなければ再利用可とする。マニフェストありでも判定の中心はケース成功状態と内部JSONのSHA。
- [`scenario_config.py:274–279`](../src/support/scenario_config.py) のケース識別子には生成器/SIFコードとHeatSolveの変更が含まれない。

**変更**

- 新方式はマニフェスト必須とし、存在しない場合は信頼しない。
- メッシュは `mesh_hash` のみならず生成器バージョン、ElmerGrid変換条件、実メッシュファイルと面接続チェックサムを検査。
- 定常再利用は実際のSIF SHA、メッシュ指紋、Solver／HeatSolve／UDF指紋、成功状態、全TESの `.state` が一致した場合のみ許可する。
- 特に **Mortarをoffへ変更した際は定常も新しい条件で解き直す**。変更前の `.result` を新しい過渡へ投入しない。
- MPIでは全ランクのrestartと全TESの状態ファイルを検査する。

### P1-03 — UDFの探索順を固定する

**根拠**: [`scenario_cli.py:28–41`](../src/support/scenario_cli.py) は新UDFを `workspace/work/udf` に置くが、[`run.py:99–122`](../src/run.py) は旧 `build/udf` も検索経路へ先に入れる。プラットフォームによって同名ライブラリの探索順が不透明になり得る。

**変更**: 新方式のUDFは絶対パスまたは明示的な専用ライブラリ領域に固定。ビルド時と実行時のABI／バイナリ指紋を確認し、旧UDFとの混在を停止条件にする。

### P1-04 — PoSTの中央パルス位置を細分化する

**根拠**: [`generate_hybrid_prism_geometry.py:697–712`](../src/generate_hybrid_prism_geometry.py) は局所メッシュをTES位置と追加 `--refine-x` に配置するが、[`scenario_cli.py:55–68`](../src/support/scenario_cli.py) はPoSTの中央パルス用の `--refine-x` を渡さない。PoSTのTES位置は x = −7, −3, +3, +7 mm、パルス中心は吸収体中央、`sigma=50 µm`、既定のグローバルサイズは50 µm。

**変更**: 物理パルスの位置をメッシュ生成へ渡し、パルス局所領域を細分化。エネルギーの離散正規化だけで満足せず、温度場と4TES電流の局所サイズ依存を調べる。位置が異なるパルスを走らせる場合は、メッシュ識別子にもその細分化設定を反映する。

### P1-05 — Stycast寸法の唯一の真値を定義する

**根拠**: [`projects/models/post_four_tes.toml`](../projects/models/post_four_tes.toml) と [単ピクセルモデル](../projects/models/single_pixel.toml) は `Stycast_dx=400[um]`、`Stycast_dy=Stycast_dx` を持つ。一方、[`scenario_config.py:230`](../src/support/scenario_config.py) は `layout.stycast_diameter` の既定値を `498[um]` としており、[`scenario_cli.py:59`](../src/support/scenario_cli.py) が後者をジオメトリ生成へ渡す。

**変更**: COMSOLの接触径と検証済み単ピクセルメッシュの物理仕様を照合し、直径を `layout` に明示するなど、幾何寸法の定義を一元化。 `Stycast_dx/dy` が本当に不要なら削除する。単純に400↔498 µmを書き換えて物理仕様まで変更しない。

### P2-01〜04 — レガシー隔離と診断性

- [`projects/templates/pipeline_template.json`](../projects/templates/pipeline_template.json) の旧 `geometries`／多数の `meshes`／`elmer_overrides` を新方式へ丸ごと持ち込まない。特に `fragment_mortar_interfaces=true` は現行ハイブリッドメッシュ生成器とは別経路であり、[旧開発記録](dual_tes_plan.md) はこのフラグの元実装をno-opだったと報告している。**現在の新方式の波形不良をこのフラグ単独に帰さない**。
- 新モデル内の `dabs_dx`、`tes_pitch`、`membrane_dz` などは依存関係を調査し、不要なら削除・非反映化。現在の [`scenario_config.py:273–278`](../src/support/scenario_config.py) は `*_dx/dy/dz` をメッシュ識別ハッシュへ広く含めるため、使われない値でメッシュ名が変わり得る。実際のジオメトリ生成引数に基づくハッシュへ変更する。
- 現行 `circuit_inner` と、レガシーの `circuit_local` / `circuit_implicit` / `circuit_parallel`、Mortar Block-Schurを同じSIFビルダー内で分岐させるのを縮小する。再現用途の旧経路は壊さない。
- 通常の `save_field=false` は容量削減のため維持できるが、`diagnostic` 実行では過渡場とTES温度、吸収体温度、界面熱流束、投入/流出熱量、4回路の反復履歴を保存できるようにする。

## 実装順序と回帰ゲート

### フェーズA：実行契約と観測性（物理条件は変えない）

1. 新方式 `show / run --dry-run` に、**有効な**Solver名、Mortar状態、時間窓、BDF次数、予測子、Step Commit、ロードするElmerSolver/HeatSolve/UDF、mesh/restart由来を表示する。
2. `[solver]` を厳密検証し、未知名は例外とする。
3. ランタイムpreflightと内容指紋を導入し、由来不明のrestart/meshを再利用しない。
4. 単ピクセル既存設定を変更せず、生成SIFと実行結果のベースラインを固定する。

**ゲートA**: 設定誤記・標準版HeatSolve誤選択・古いDLL・不明なrestartを事前に検出できる。既存の単ピクセル基準ケースに意図しない変更がない。

### フェーズB：完全共形ルートへ統一（物理影響を分離）

1. PoSTの各界面で、全節点共有／境界と体積要素のトポロジー／熱流束を先に検査。
2. 定常と過渡の `apply_mortar_bcs=False` を明示。新しい識別子で**定常から**再計算。
3. 旧ケースと、同一メッシュで「Mortarあり/なし」を比較して界面拘束による差を測定。
4. 単ピクセルのSIF/波形とCOMSOLの比較結果を再確認する。

**ゲートB**: Mortarなしで全界面が連続し、漏れた熱源や絶縁された部品がなく、単ピクセル基準を破壊していない。旧非共形Mortarケースは従来どおり生成可能。

### フェーズC：PoST過渡波形の資格確認

1. PoST 2案件を75 ms観測へ変更し、BDF2とStep Commitの採用方針を固定。
2. パルスなしで定常restart後の4回路電流・温度ドリフトを確認。
3. パルスありの4回路電流・TES温度・吸収体温度・熱収支を保存。
4. パルス位置のメッシュ細分化を変え、波形の空間収束を評価。
5. 同じメッシュ・初期状態でMUMPS対HYPREを比較し、最後にCOMSOLを対照する。

**ゲートC**: 全4チャンネルの波形がパルス前から減衰まで再現可能。数値差と物理差が独立に説明できる。COMSOLとの許容誤差は基準データと観測時間を明記して設定し、比較前に合否閾値を固定する。

## 必須自動テスト（追加・拡張）

主に [`tests/test_scenario_config.py`](../tests/test_scenario_config.py)、[`tests/test_run_restart_reuse.py`](../tests/test_run_restart_reuse.py) とSIF/メッシュ専用テストで実施。

- [ ] 現行の単ピクセル／PoSTの定常・過渡SIFにMortar/Projectorが出現しない。
- [ ] 4つのStycast–Pb接触を含む全熱接続が実メッシュで共有節点／接続性の条件を満たす。
- [ ] レガシー非共形ケースはMortarありで従来のSIFを生成できる。
- [ ] PoST `window` 未指定時に意図しない短時間・BDF1へサイレントに落ちない。
- [ ] BDF次数、予測子、Step Commitが `show` とSIFで一致し、選定した回路更新がテストで固定される。
- [ ] `[solver] linear_system="mumsp"` や未知のキーを拒否する。
- [ ] 標準版HeatSolve、古いUDF、互換性のないツールチェーンをpreflightで拒否する。
- [ ] SIF生成コード／実行バイナリ／mesh変更で定常restart再利用を無効化する。履歴不明のマニフェストも拒否する。
- [ ] 不使用の幾何パラメータの変更ではメッシュ内容と識別子が変わらず、実使用寸法の変更では再生成する。
- [ ] PoST中央熱源の離散積分、局所分割サイズ収束、4TES電流・温度のチャンネル対応を検証する。
- [ ] パルスなしのrestartで、4チャンネルの不自然なベースライン緩和が起きない。
- [ ] 一連の結果に完全な実行環境・メッシュ・初期状態・Solver設定の指紋が残る。

## 触らずに維持すべきもの

- COMSOLとの照合に使った**物性値・TES回路値・実際の接触寸法**。設定の掃除と物理変更を混同しない。
- レガシー専用コマンド（`python main.py legacy ...`）と既存の旧非共形ケースの再現能力。
- 意図的に使われる `circuit_inner`、非線形連成、`TES Parallel Power`、単ピクセルの検証済みh8メッシュ。
- 一時的な診断設定（例えばMortarのA/B検証）は本番既定値へ戻さず、明示的な実験ケースとして保存する。

## 実行結果についての留保

本計画に記載した「生成され得る」「不整合」「互換性リスク」は、参照コミットの静的コード解析から判断した。**この文書自体はPoSTの波形不良の根本原因を実証するものではない**。原因を断定するには、生成済みSIF、実メッシュとロードされたHeatSolve、定常restart、4TESの時系列、COMSOL参照データを同条件で突き合わせる必要がある。
