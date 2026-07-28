# 吸収体テトラ／積層部プリズムのハイブリッドメッシュ：現状

更新日: 2026-07-25

## 目的

MPI計算での要素数を抑えつつ、吸収体中央に入射するパルスを十分に扱えるよう、以下の別メッシュを試作した。

- 吸収体 `abs`: 四面体要素
- TES、Stycast、膜、SiO2、Si、SiNx: 三角形断面を厚み方向に押し出したプリズム要素

既存の全四面体メッシュ `mesh_refined_3x` とそのMPI分割メッシュは変更していない。

## 既存メッシュとの比較

| 項目 | 既存 `mesh_refined_3x` | 新規ハイブリッド |
|---|---:|---:|
| 体積要素数 | 312,997（全て四面体） | 127,037 |
| 四面体 | 312,997 | 42,803（吸収体） |
| プリズム | 0 | 84,234（積層部） |
| 節点数 | 63,172 | 60,093 |

体積要素数は約59%削減できた。薄膜の厚みはプリズム押し出しで表現するため、1--20 umの膜が吸収体全体の三次元テトラ細分化を要求しない。

## 生成物

- 生成器: `generate_hybrid_prism_geometry.py`
- Elmerメッシュ: `mesh_hybrid_abs_tet_layers_prism_conformal/`
- 比較用プロジェクト: `elmer_project_hybrid_prism.json`
- ケース定義生成器: `scripts/prep/prepare_hybrid_prism_case.py`

プロジェクト生成:

```powershell
python scripts/prep/prepare_hybrid_prism_case.py
python run.py case_tes_pulse_hybrid_prism_fast_compare --project elmer_project_hybrid_prism.json --dry-run
```

メッシュ生成・変換:

```powershell
python generate_hybrid_prism_geometry.py elmer_project_comsol_timegrid.json
& 'C:\Program Files\Elmer 26.1-Release\bin\ElmerGrid.exe' 14 2 `
  gmsh/project_hybrid_prism.msh -merge 1e-10 `
  -out mesh_hybrid_abs_tet_layers_prism_conformal
```

ElmerGrid後の `mesh.header` は、体積要素として `504`（四面体）と `706`（プリズム）を含む。

## 熱源設定

現在のパルス熱源は、吸収体の体積重心を中心とする三次元ガウス分布である。

- 中心: `(x, y, z) = (0, 1.000 mm, 0.56216 mm)`
- 空間幅: `sigma = 50 um`
- パルス: 20.02 ms開始、1 nsの矩形時間窓
- エネルギー: 1332 keV

ハイブリッドメッシュでも吸収体重心は同じである。パルスを半径50 umの一様球に変更することは可能だが、その場合は吸収体中心近傍の局所テトラ細分化（目安10--20 um）が必要になる。

## 検証済み

標準Elmer 26.1で、ハイブリッドメッシュの定常ケースは完走した。

- ケース: `case_tes_steady_hybrid_prism`
- 結果: `mesh_hybrid_abs_tet_layers_prism_conformal/case_tes_steady_hybrid_prism.result`
- 回路状態: `mesh_hybrid_abs_tet_layers_prism_conformal/case_tes_steady_hybrid_prism.state`

標準Elmerでの完走により、Gmsh生成、ElmerGrid変換、プリズム／テトラ混在要素、材料・モルタル境界の基本的な読み込みは確認できている。

## 解決済み: カスタム実行時のUDF ABI不整合

初期に使用していた `tools/elmer-hypre/install` は、`libelmersolver.dll` が旧 `build`、`HeatSolve.dll` が `build-gfortran16` 由来で混在していた。さらにプロジェクト直下の旧 `tes_transient_heat_source_t0.dll` を読み込んだ場合、同一SIFは `HeatSolve: Assembly` 直後でSIGSEGVとなった。

gfortran 16 / UCRT64で、Elmer本体・`libelmersolver`・`HeatSolve`を新規の同一ツリーへReleaseビルドし、同じinstallの `elmerf90` でUDFも再ビルドしたところ、ハイブリッド定常ケースは完走した。したがって、プリズム混在要素またはモルタル境界そのものの非互換性ではなく、旧UDF DLLとのABI不整合が原因である。

- matching Release build/install: `build-phase1-gfortran16-release-retry5` / `install-phase1-gfortran16-release-retry5`
- matching UDF SHA256: `630303a0b9d8d326392ba17e7c28d7ce607e06cfea7fff5948a4474cdc5083b4`
- matching-UDF定常: exit 0、188 s
- 標準Elmerとの差: 温度の最大絶対差 `1.95e-7 K`、TES state温度差 `2.57e-7 K`
- build/provenance: `artifacts/hybrid_prism_diagnostics/build-release-20260724-201700-retry5/`
- matching-UDF実行・比較: `artifacts/hybrid_prism_diagnostics/20260724-202037/`

再現例（標準結果を再実行せずcustomだけを確認）:

```powershell
.\scripts\support\run_hybrid_prism_ab.ps1 -CustomOnly `
  -CustomSolver ..\tools\elmer-hypre\install-phase1-gfortran16-release-retry5\bin\ElmerSolver.exe `
  -CustomRuntimeBin C:\msys64\ucrt64\bin `
  -UdfDll artifacts\hybrid_prism_diagnostics\udf-matching-release-retry5-20260724-201900\tes_transient_heat_source_t0.dll `
  -TimeoutSeconds 300
```

## MPI回帰の現状

`repart_x` 4-rank partitionは、体積要素数と `504` / `706` の要素型別個数を保存しており、定常・短時間過渡ともexit 0で完走する。serial-tight定常場を各rankへ写像した値は、64,759個すべてで元のserial値との差が `0 K` である。stateファイルもSHA256一致を確認した。

Phase13では `TESInnerCircuitUpdate` に opt-in の `"TES Inner Circuit Step Commit"` を追加した。各時間ステップの境界で、直前の最終sweep温度、前ステップでcommit済みの電流、前ステップ幅からBackward Euler式を再評価し、その確定値を次ステップの履歴状態に使う。従来の「最後に実行された非線形sweepの試行値をそのまま履歴へ渡す」経路と、sweep回数に応じたmidpoint補正をstep-commit経路から除いた。新ステップ内の試行電流・Joule powerは従来どおり熱方程式と反復する。対象ケースは自然収束であり、`Nonlinear System Min Iterations` は設定していない。

共通20 ms restartを使うstep-commitの1-rank／4-rank比較は次のとおり。

| 指標 | 差 | 基準 | 判定 |
|---|---:|---:|:---:|
| 基線電流 | 0.0332% | 1% | PASS |
| パルス波高 | 1.7347% | 2% | PASS |
| ピーク時刻 | 0 s | 10 us | PASS |
| baseline補正後の最大波形差 | 1.8226% | 2% | PASS |
| raw最大電流差 | 2.4809% | 2% | FAIL |

独立した4基準（基線・波高・時刻・baseline補正波形）はPASSである。一方、baseline offsetも含む診断値であるraw最大電流差は2%に未達のため、rawを含む5基準の全面PASSとはしない。反復回数は47ステップ中2ステップで異なる（1-rank: 1回×3、2回×44、4-rank: 1回×2、2回×44、3回×1）。step境界の履歴値は決定的に再評価されるが、ステップ内の熱・回路連成反復数は完全には一致していない。

性能はstep-commit自然収束で1-rank 236.52 s／4-rank 144.29 s、並列加速は1.64xである。min3の1-rank 366 s／4-rank 244 s（1.50x）より速く、人工的な最小反復数を要求しない。従って、以後のMPI回帰の推奨設定はstep-commit自然収束とする。`-partnobcoptim` partitionはMPI組立でSIGSEGVとなったため不採用である。Phase12のmin5もserial 560.76 s／MPI 377.7 s、corrected差3.243%、raw差3.912%へ悪化したため不採用である。

実装・再現用の固定成果物は次のとおり。

- 実装: `tools/elmer-hypre/src/fem/src/modules/HeatSolve.F90`、`scripts/support/build_cases.py`、`scripts/prep/prepare_hybrid_prism_case.py`
- matching build/install: `tools/elmer-hypre/build-phase13-step-commit` / `tools/elmer-hypre/install-phase13-step-commit`
- matching UDF: `artifacts/hybrid_prism_diagnostics/phase13_udf/tes_transient_heat_source_t0.dll`
- 結果・来歴: `artifacts/hybrid_prism_diagnostics/phase13/comparison.json`、`iteration_distribution.json`、`run_history.json`

主要な検証結果:

- `artifacts/hybrid_prism_diagnostics/phase10/mapping_validation_v2.json`
- `artifacts/hybrid_prism_diagnostics/phase11_comparison.json`
- `artifacts/hybrid_prism_diagnostics/phase11_iteration_distribution.json`
- `artifacts/hybrid_prism_diagnostics/phase12_comparison.json`
- `artifacts/hybrid_prism_diagnostics/phase12_iteration_distribution.json`
- `artifacts/hybrid_prism_diagnostics/phase13/comparison.json`
- `artifacts/hybrid_prism_diagnostics/phase13/iteration_distribution.json`
- `artifacts/hybrid_prism_diagnostics/phase13/run_history.json`
- `artifacts/hybrid_prism_diagnostics/phase14/all_tet_vs_hybrid_step_commit.json`
- `artifacts/hybrid_prism_diagnostics/phase14/body_integrals_step_commit_tight.json`
- `artifacts/hybrid_prism_diagnostics/phase14/all_tet_iteration_summary.json`
- `artifacts/hybrid_prism_diagnostics/phase14/hybrid_iteration_summary.json`
- `artifacts/hybrid_prism_diagnostics/phase14/runtime_summary.json`
- `artifacts/hybrid_prism_diagnostics/phase15/mesh_validation.json`
- `artifacts/hybrid_prism_diagnostics/phase15/body_integrals_stack25_vs_all_tet.json`
- `artifacts/hybrid_prism_diagnostics/phase15/steady_comparison.json`
- `artifacts/hybrid_prism_diagnostics/phase15/all_tet_vs_stack25_step_commit.json`
- `artifacts/hybrid_prism_diagnostics/phase15/stack25_iteration_summary.json`
- `artifacts/hybrid_prism_diagnostics/phase15/pulse_improvement.json`
- `artifacts/hybrid_prism_diagnostics/phase15/runtime_summary.json`
- `artifacts/hybrid_prism_diagnostics/phase16/mesh_validation.json`
- `artifacts/hybrid_prism_diagnostics/phase16/body_integrals_stack17_vs_all_tet.json`
- `artifacts/hybrid_prism_diagnostics/phase16/steady_comparison.json`
- `artifacts/hybrid_prism_diagnostics/phase16/all_tet_vs_stack17_step_commit.json`
- `artifacts/hybrid_prism_diagnostics/phase16/stack17_iteration_summary.json`
- `artifacts/hybrid_prism_diagnostics/phase16/pulse_improvement.json`
- `artifacts/hybrid_prism_diagnostics/phase16/runtime_summary.json`
- `artifacts/hybrid_prism_diagnostics/phase17/mesh_validation.json`
- `artifacts/hybrid_prism_diagnostics/phase17/radial_refinement_diagnostic.json`
- `artifacts/hybrid_prism_diagnostics/phase18/mesh_validation.json`
- `artifacts/hybrid_prism_diagnostics/phase19/mesh_validation.json`
- `artifacts/hybrid_prism_diagnostics/phase19/all_tet_vs_phase19_step_commit.json`
- `artifacts/hybrid_prism_diagnostics/phase19/pulse_iteration_summary.json`
- `artifacts/hybrid_prism_diagnostics/phase19/body_integrals_phase19_vs_all_tet.json`
- `artifacts/hybrid_prism_diagnostics/phase19/final_summary.json`

## 全四面体基準との差

Phase14では、Phase13と同じmatching runtime/UDF、step-commit回路、20 ms restart、8段時間格子を用いた。全四面体側は新規のtight定常解から開始しており、初期化経路と回路モードを揃えた比較である。全四面体定常は353 s、過渡は922 s、ハイブリッド過渡は236.52 sであり、ハイブリッドは約3.90倍高速である。

| 指標（全四面体を基準） | 差 | 基準 | 判定 |
|---|---:|---:|:---:|
| 基線電流 | 11.0105% | 1% | FAIL |
| パルス波高 | 6.1801% | 2% | FAIL |
| ピーク時刻 | 100 us | 10 us | FAIL |
| baseline補正後の最大波形差 | 15.4224% | 2% | FAIL |
| raw最大電流差 | 208.8361% | 2% | FAIL |

基線電流は全四面体148.153780 uA、ハイブリッド164.466321 uA、パルス波高は順に7.814483 uA、8.297423 uAである。従って、現時点で全四面体／ハイブリッド間のメッシュ収束は達成していない。

定常体積積分では、TES平均温度は全四面体0.168438210 K、ハイブリッド0.167983313 Kで、全四面体が+0.454897 mK高い。吸収体とStycastもほぼ同じ+0.455 mKの動作点シフトを示す。TES--Stycast平均温度差は全四面体-0.351687 uK、ハイブリッド-0.399102 uK、TES--Membrane_SiNx差は5.234690 mK、4.693234 mKである。Stycast体積だけは+0.144229%異なり、他bodyの体積差は数値丸めの範囲である。従って、均一な約0.455 mKの動作点移動と約10.3%のTES--膜勾配不足は、bulk体積の欠落ではなく、積層部／膜の面内離散化または接触の実効熱抵抗を主因として示唆する。

過渡の反復分布は全四面体が `{1回: 2, 2回: 44, 8回: 1（step 42）}`、ハイブリッドが `{1回: 3, 2回: 44}` である。なおPhase14で、HeatSolveの`Element%BodyId`は物理target IDではなくSIFのBody ordinalであることを確認した。ケース生成器は全四面体TESをordinal 2、ハイブリッドTESをordinal 8として出力するよう修正し、非連続physical target IDを使う回帰テストを追加済みである。

Phase15では中央stack／接触柱だけを、中心からx/y各±0.4 mm、25 umのBox fieldで細分化し、50 umの遷移を設けた。Boxは吸収体底面の接触面まで届くため、その界面直近のtetraも細分化されるが、吸収体bulkは50 umを維持する。stack25 meshは67,069ノード、147,995体積要素（tetra 58,558、prism 89,437）であり、legacy hybridより+16.50%、全四面体より-52.72%の体積要素数である。局所disk再テッセレーションによりStycast体積はlegacy比+0.478689%、全四面体比+0.623608%となるため、この体積差は幾何tessellationの交絡要因として扱う。

stack25定常は159.45 sで完走し、最終TES stateは0.168325831 K、152.122818 uAである。全四面体に対するlegacy→stack25のgap閉鎖率は、TES平均温度75.30%、state電流75.67%、TES--Membrane_SiNx勾配62.74%である。中央stackの解像度が動作点差の主要因であることを強く示すが、これだけで収束達成とはしない。

stack25のstep-commit過渡は381.40 sで完走した。

| 指標（全四面体を基準） | 差 | 基準 | 判定 |
|---|---:|---:|:---:|
| 基線電流 | 2.6789% | 1% | FAIL |
| パルス波高 | 1.5242% | 2% | PASS |
| ピーク時刻 | 0 s | 10 us | PASS |
| baseline補正後の最大波形差 | 7.2400% | 2% | FAIL |
| raw最大電流差 | 50.8896% | 2% | FAIL |

Phase14 legacy hybridに比べ、基線・波高・raw差は約75.7%／75.3%／75.6%、baseline補正後差は53.1%、ピーク時刻差は100%低減した。反復分布は `{1回: 2, 2回: 34, 3回: 10, 8回: 1（step 42）}` である。したがって、stack25は動作点とピーク位置を大幅に改善した一方、baseline補正波形は未収束であり、吸収体パルス中心の局所細分化も必要と考えられる。

Phase16では同じ中央stack／接触柱を16.7 um（遷移50 um）へ細分化した。stack17 meshは76,630ノード、170,960体積要素（tetra 73,282、prism 97,678）であり、stack25比+15.52%、全四面体比-45.38%である。局所disk再テッセレーション後のStycast体積はstack25比+0.091399%、全四面体比+0.715577%であり、引き続きtessellation由来の交絡要因として記録する。

stack17定常は184.84 sで完走し、最終TES stateは0.168481471 K、146.628193 uAである。全四面体に対して温度は+43.452 uK、電流は-1.526379 uAへ僅かにovershootした。一方、TES--Membrane_SiNx平均温度差は5.215751 mK（全四面体5.234690 mK）で、legacyから全四面体までのgapを96.50%閉じた。TES--Stycast平均温度差は-0.338218 uK（全四面体-0.351687 uK）、Stycast--absorber差は-0.000053 uK（全四面体-0.000148 uK）である。したがって、中央stackの定常動作点はさらに近づいたが、全てのbody／界面指標で収束を宣言する段階ではない。

stack17のstep-commit過渡は721.43 sで完走した。

| 指標（全四面体を基準） | 差 | 基準 | 判定 |
|---|---:|---:|:---:|
| 基線電流 | 1.0298% | 1% | FAIL |
| パルス波高 | 0.1990% | 2% | PASS |
| ピーク時刻 | 100 us | 10 us | FAIL |
| baseline補正後の最大波形差 | 6.3849% | 2% | FAIL |
| raw最大電流差 | 25.9091% | 2% | FAIL |

stack25に対して、stack17は基線・波高・baseline補正後差・raw差をそれぞれ1.0298%・0.1990%・6.3849%・25.9091%まで下げたが、ピーク時刻は再び100 usずれた。反復分布は `{1回: 2, 2回: 26, 3回: 17, 10回: 1（step 43）, 12回: 1（step 42）}` であり、stack25の最大8回より剛性が増している。要素数はstack25比+15.52%だが、過渡時間は+89.15%（381.40 sから721.43 s）である。

次の吸収体パルス中心の局所細分化にはstack17を基準メッシュとして採用する。これはstack25よりピーク時刻だけは悪化したものの、基線は1%閾値直前まで縮小し、波高・補正波形・raw波形の全てで全四面体に近づいたためである。ただし全基準PASSではなく、反復剛性と計算時間の増加も大きい。従って、これ以上の中央stack全体の細分化は行わず、吸収体中心だけを局所的に細分化して残る波形差とピーク時刻を評価する。

Phase17ではstack17に吸収体中心16.7 um／半径150 umのBall fieldを重ねたが、943,968体積要素（tetra 846,290）となり250k上限を大幅に超えた。半径350 um外（吸収体体積の約74%）の等価tetra edge長中央値も33.93 umから17.91 umへ低下しており、局所細分化ではなく全域への波及であるため不採用とした。Phase18では`Mesh.MeshSizeExtendFromBoundary=0`、25 um／半径75 umへ縮小したが、なお421,923要素で上限超過、350 um外の中央値も23.92 umであった。この二つはsolverを実行していない。

Phase19では同じboundary-extension無効化の下で35 um／半径50 umとし、223,718体積要素（tetra 125,972、prism 97,746）へ収めた。これは全四面体比-28.52%、stack17比+30.86%である。350 um外には依然74.35%の吸収体体積があり、等価edge長中央値30.77 um（stack17 33.93 um）と軽微な全域波及は残る。physical nameとTES Body ordinal 8を保持し、body体積はStycastの全四面体比+0.715577%を除き丸め範囲である。

Phase19定常は191.06 s、最終TES状態は0.168460630 K、147.359511 uA、15.023108 mOhm、3.262242e-10 Wである。全四面体に対して温度は+22.612 uK、電流は-0.795061 uAとなった。TES--Membrane_SiNx平均温度差は5.188952 mK（全四面体5.234690 mK）で、legacyからのgapを91.55%閉じた。state pathは当初137文字でUDFの`CHARACTER LEN=128`を超え、122-byteのtruncated filenameとなった。既存stateを同じmesh directory内の`phase19_steady.state`へ移動し、single/dualとも128文字超をbuilderで拒否する検証を追加して修正した。定常は再実行していない。

Phase19 step-commit過渡は801.96 sで完走した（全四面体922 s比-13.02%）。

| 指標（全四面体を基準） | 差 | 基準 | 判定 |
|---|---:|---:|:---:|
| 基線電流 | 0.5371% | 1% | PASS |
| パルス波高 | 0.0244% | 2% | PASS |
| ピーク時刻 | 100 us | 10 us | FAIL |
| baseline補正後の最大波形差 | 6.5114% | 2% | FAIL |
| raw最大電流差 | 16.6933% | 2% | FAIL |

基線電流は全四面体148.153780 uA、Phase19は147.358117 uA、波高は7.814483 uA、7.812576 uAである。反復分布は`{1回: 2, 2回: 24, 3回: 19, 10回: 1（step 43）, 12回: 1（step 42）}`であり、stack17と同じ最大反復数である。Phase19は基線と波高を初めて同時に閾値内へ入れ、raw差もstack17の25.9091%から16.6933%へ下げた。一方でピーク時刻、baseline補正後波形、raw波形の3基準は未達であり、メッシュ収束の宣言はしない。

したがって、Phase19を次のCOMSOL波形比較のbounded hybrid candidateとする。全四面体より少ない要素数と短い過渡時間で基線・波高を満たすためである。これ以上の吸収体global refinementは行わず、COMSOLと同じ時刻・基線補正・パルス位置で、ピーク100 us差と波形残差をまず切り分ける。残差がCOMSOL比較でも持続する場合に限り、熱源時間位置、absorber中心の局所field遷移、Stycast／膜の面内分割・接触熱抵抗を一因子ずつ評価する。

## 次の作業

1. Phase19を固定し、COMSOLと同一の時間格子・基線補正・パルス中心で波形を比較する。特に100 usのピーク時刻差と6.5114%の補正波形差を独立に確認する。
2. COMSOL比較でも残差が持続する場合に限り、absorber中心fieldの遷移、Stycast／膜の面内分割、接触熱抵抗を一因子ずつ評価する。250k要素上限を維持し、global absorber refinementは行わない。
3. 5基準すべてが閾値内になるまでは、全四面体／ハイブリッド間のメッシュ収束を宣言しない。
