# TES モデルの拡張ガイド

> **Legacy workflow:** この文書は `cases/*.toml` + `tes_sim.py` の旧拡張経路を説明しています。通常CLIからは `python main.py legacy toml ...` を入口にしてください。新規計算は `projects/*.toml` + `projects/models/*.toml` を使い、多段 override は行いません。新方式は [シミュレーションの始め方](getting_started.md) と [project / model TOML 設計](toml_unification_plan.md) を参照してください。

単ピクセル以外の構成（2TES、正方形の4隅・三角形の3隅に TES を置く構成など）や、材料・配線・回路の変更を行う実装者向けの手引きです。高速化済みの計算経路（Phase24 HYPRE、6ランク MPI）をそのまま使えるように、どこを書けば何が変わるかをまとめています。

## 1. 全体の流れ

```
メッシュ（Gmsh など、形は自由）
   │  python tes_sim.py mesh-import foo.msh --name mesh_foo [--roles roles.toml]
   ▼
メッシュ登録（artifacts/phase24_conformal_hybrid/project.json の "meshes"）＋役割表
   │  cases/<名前>.toml（メッシュ名・役割表・TES ごとのパラメータ・パルス・時間窓）
   ▼
python tes_sim.py run cases/<名前>.toml
   │  分割 → 定常 → 過渡 → 要約（runs/<名前>/summary.json）
   ▼
scripts/support/build_cases.py が SIF を生成 → HeatSolve（N 回路）で計算 → results/<case>/
```

- 形状はメッシュが決めます。コード側は、本体と境界の **名前** と **役割表** しか見ません。
- 役割表を書かなければ、従来の命名規則（`TES`、`abs`、`SiO2_2`… と接尾辞 `_L`／`_R`／`_T<k>`）で役割を推測します。

## 2. 役割表（roles）

メッシュ登録（`project.json` の `meshes.<名前>.roles`）か、ケースファイルの `[roles]` に書きます。ケースファイルの記述が優先されます。

```toml
[roles]
tes_circuits    = ["TES_left", "TES_right"]   # TES 回路として扱う本体。並び順が回路番号 1..N
pulse_bodies    = ["abs"]                     # パルスを与える本体（パルス中心の自動計算と正規化にも使う）
bath_boundaries = ["bath"]                    # 浴温度（T_bath）に固定する境界

[roles.bodies]                                # 本体名（fnmatch のワイルドカード可、上から最初に一致したもの）→ 材料
"TES_*"      = "TES"
"Membrane_*" = "Membrane"
"Si*"        = "Si"
"wire_*"     = "Al"                            # project.json の "materials" に追加した材料
"abs"        = "Pb"
```

- 材料名は `project.json` の `"materials"` のキーです。標準の7材料（Pb, TES, Stycast, SiO2, Si, SiNx, Membrane）は SIF の Material 1〜7 に固定です。追加した材料は 8 番以降に並びます。
- 役割表で書かなかった項目は命名規則で補います。
- 名前がメッシュにない場合は、SIF を作る前にエラーで止まります（`validate_roles`）。

## 3. よくある拡張の手順

### 3.1 TES の配置を変える（2TES、4隅、3隅など）

1. Gmsh などでメッシュを作ります。条件は次の3つです。
   - TES 回路にする本体を、回路ごとに1つの物理グループ（本体）にする。
   - 浴に固定する面を境界の物理グループにする。
   - 要素は四面体（504）とプリズム（706）にする。これ以外の要素は汎用の遅い経路で組み立てられます。
2. 界面は節点を共有させる（共形）のが基本です。共形にしない場合は、mortar 境界条件の設定が別途必要です。参考として `generate_hybrid_prism_geometry.py` の `--conformal-tes-stack`／`--conformal-abs` の作り方（TES の輪郭と Stycast 円盤を各層に刻み、吸収体側は周期写像で節点を一致させる）があります。直線上に並べるだけなら `--tes-x-positions` で作れます。
3. `python tes_sim.py mesh-import foo.msh --name mesh_foo --roles roles.toml` で登録します。
4. ケースファイルを書き、`python tes_sim.py run` で実行します。

TES の数は最大8です。発熱の UDF の入口 `TESParallelHeatSource1..8`（`tes_parallel_circuit.f90`）で決まっていて、関数をコピーすれば増やせます。

### 3.2 TES ごとにパラメータを変える

ケースファイルに、回路の数だけ `[[tes]]` を並べます（順番は `tes_circuits` と同じ）。書いたキーだけが上書きされます。

```toml
[[tes]]            # 回路 1
I_bias = "700[uA]"
[[tes]]            # 回路 2
alpha = 240
R_0 = "16[mOhm]"
```

使えるキー（`build_cases.TES_CIRCUIT_KEYS`）: `I_bias, R_sh, L_tes, R_0, R_min, alpha, beta, I_0, T_c, T_0, TES_volume`。値は SI の数値か、単位付きの式です。

### 3.3 温度依存性を別の場所に付ける

`project.json` の `"materials"` で、対象の材料の `k`（`rho`、`cp` も同様）の `expression` に `T` を書くだけです。

```json
"Stycast": {"k": {"expression": "2.69e-6*(T/0.17)**2"}, ...}
```

- `T` を含む式は、3次スプラインの温度表（`Real Cubic`、既定 0.05〜1 K を対数等間隔に400点）として SIF に出力されます。UDF は不要です。
- 範囲と点数は、ケースファイルの `materials_temperature_tables = {Stycast = {t_range = [0.1, 0.3], points = 400}}` で変えられます。
- メンブレンだけは、検証済みの UDF（`tes_membrane_conductivity.f90`、べき乗則）を使います。
- **熱伝導率 k(T)** は、高速組み立て（プリズム）で何材料でも扱えます。節点の k は「材料×節点」ごとに保持しています（`HeatSolve.F90` の `Phase24DynK*`、`Phase24EvaluateDynK`）。非線形反復の終了判定（熱伝導率の遅れ）と集中化前処理も、すべての温度依存材料を対象にしています。
- **比熱 cp(T)・密度 ρ(T)** は正しく計算されますが、その材料の要素は汎用の遅い経路になります。高速化するには、`Phase24AssemblePrismElements` の質量（`Phase24PrismRhoCp`、`Phase24PrismMass`）を、動的な k と同じ要領で節点値から毎回作り直す処理を追加します。

### 3.4 アルミ配線を置く

1. メッシュに配線の本体を加えます。膜の上の薄膜なら、`generate_hybrid_prism_geometry.py` の `_imprint_tes_footprint` と同じ要領で配線の帯の輪郭を層に刻み、別の本体として押し出すと、節点を共有したまま入れられます。
2. `project.json` の `"materials"` に `Al` を追加します。超伝導 Al の熱伝導率は温度依存が強いので、`T` の式で書きます（3.3 のとおり表になります）。
3. 役割表の `bodies` で、配線本体を `Al` に割り当てます。

この方法で扱えるのは熱の経路としての配線だけです。配線の電気抵抗や、配線に流れる電流による発熱は含まれません（3.5 を参照）。

### 3.5 TES 内を実際に流れる電流を計算する（難易度高）

現在の回路は集中定数（TES 1個に電流1つ）です。電流は `TESInnerCircuitUpdateOne` で回路方程式（シャント抵抗、インダクタ、R(T, I)）から求め、発熱 P = I²R を TES 本体の体積に一様に与えています。

TES 内の電流分布を解くには、次の作業が必要です。

1. TES（と配線）本体で電位の方程式 ∇·(σ(T) ∇φ) = 0 を解くソルバーを追加する（Elmer の `StatCurrentSolve` が候補）。電極となる境界の定義も要る。
2. 局所の電気伝導率 σ(T, j) を、転移のモデル（R(T, I) の局所版）から与える。
3. ジュール発熱 σ|∇φ|² を、熱方程式の体積発熱として毎反復で受け渡す（今の `TESParallelHeatSource` の置き換え）。
4. 外部回路（バイアス電流、シャント、インダクタ）は、電極間の電圧と総電流の関係として境界条件に組み込む（今の回路方程式の拡張）。
5. 熱と電気の連成を非線形反復に入れ、終了判定（`Phase24 Coupled Temperature Tolerance`）を電気側の残差にも広げる。

入口になるのは、`HeatSolve.F90` の `TESInnerCircuitUpdate`（反復ごとの呼び出し位置）と、`build_cases.py` の `body_force_blocks`（発熱の与え方）です。

## 4. 数値設定の目安

| 項目 | 目安 | 理由 |
|---|---|---|
| プリセット | `fast`（既定） | 75 ms を単ピクセルで約5分。精度は COMSOL との差 0.09 / RMSE 0.03 µA で、厳しい設定と同等 |
| 線形許容値 | TES 1個は 3e-11、複数は 1e-11（`fast` が自動で選択） | 系が大きいほど ‖b‖ に対する相対誤差が TES 温度の揺らぎとして効く。非線形の判定（1e-7 K）を下回らないと、反復が2状態を往復して上限まで回る |
| 定常計算 | 回路ニュートン＋ `Phase24 Coupled Temperature Tolerance`（`fast` の定常に含まれる） | 従来の緩和法（1回で残差の4%だけ修正）は温度場の変化が小さいだけで止まり、線形許容値次第で定常電流が約0.07 µA ずれた。過渡計算はこの定常解から再開するので、そのずれがパルス前の電流とピークにそのまま残る。新設定では線形許容値によるずれは0.006 µA 以内で、計算も約3倍速い |
| MPI ランク | 6（物理コア数） | 12（SMT）は最初の解きから進まなかった。分割は `tes_sim.py` が自動で作る |
| メッシュ | TES・メンブレン周りの局所サイズ 8 µm 程度 | h8 で熱コンダクタンス G が収束域 |

## 5. 変更したときの確認手順

1. **単ピクセルが変わっていないこと**
   - `tes_sim.py show cases/single_pixel_h8.toml` で設定を確認する。
   - SIF 生成を dry-run し、`artifacts/phase24_conformal_hybrid/transient_75ms/mpi6_best_fag1.sif` と比べる（同一であること）。
2. **メッシュ生成を変えたとき**: 既存の h8 メッシュを作り直し、`gmsh/project_hybrid_fullconf_h8.msh` とバイト比較する。
3. **計算結果**: パルス後 1 ms の区間（`1ms_bdf2h15`）で旧版と比べる。6ランクの計算は実行ごとに最大約 0.002 µA 揺らぐので、それ以下なら同じと見なす。

## 6. コードの所在

| 役割 | 場所 |
|---|---|
| CLI、プリセット | `tes_sim.py` |
| ケースファイルの例 | `cases/*.toml` |
| SIF 生成（役割表、材料、回路定数） | `scripts/support/build_cases.py`（`set_active_roles`、`resolve_bodies`、`material_property_lines`、`numbered_inner_constants_block`） |
| パルスの中心と正規化 | `scripts/support/mesh_quantities.py`（`set_pulse_bodies`） |
| 回路（N 個）、収束判定、刻み幅制御 | `tools/elmer-hypre/src/fem/src/modules/HeatSolve.F90`（`TESCircuitState_t`、`TESInnerCircuitUpdateOne`、`TESKey`） |
| 高速組み立て（温度依存の k） | 同上（`Phase24BuildMetadata`、`Phase24EvaluateDynK`、`Phase24AssemblePrismElements`） |
| TES の発熱 | `tes_parallel_circuit.f90`（`TESParallelHeatSource`、`TESParallelHeatSource1..8`） |
| 線形ソルバー（HYPRE、前処理の使い回し） | `tools/elmer-hypre/src/fem/src/SolveHypre.c`、`SParIterSolver.F90` |
| 実行（WSL、MPI） | `scripts/support/run_phase24_gpu_elmergpu.py`、`run.py` |

Elmer 本体の変更は `artifacts/phase24_conformal_hybrid/elmer_native_source.patch`（上流 9916c3f05 からの累積）にまとめてあります。WSL のビルド手順は `scripts/support/build_elmer_hypre_cpu_elmergpu.sh` です。
