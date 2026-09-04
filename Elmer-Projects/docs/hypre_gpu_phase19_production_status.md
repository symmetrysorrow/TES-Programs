# Phase 19 HYPRE GPU / mortar production status

更新日: 2026-09-04

## 結論

Phase A のクリーン適用と CPU HYPRE ビルドは完了した。Phase B の 4 ベクトル診断も、以前の 1 ベクトルだけの実行から前進し、全段階が有限値になった。

ただし、独立 SciPy/SuperLU oracle との既存の厳格ゲート（vector 相対誤差 `<=1e-14`、action 相対誤差 `<=1e-10`）は未達である。したがって指示どおり、CPU lower/full、GPU、短時間 transient、MPI へは進めていない。反復回数上限到達や solver residual のみを合格とはしていない。

## Phase A: clean base / build

- source base: `5a8de867068be0568f09af40fb90fee300dfbede`
- feature patch: [hypre_gpu_phase19_schur_feature.patch](hypre_gpu_phase19_schur_feature.patch)
- `git apply --check`: clean
- build: isolated source worktree + Ninja, MPI/HYPRE/MUMPS/UMFPACK enabled
- HYPRE library: `tools/hypre-cuda-install/lib/libHYPRE.so` (CPU execution mode)
- executable: `tools/elmer-phase19-feature-cpu-install/bin/ElmerSolver_mpi`
- optional serial SuperLU diagnostic wrapper is enabled by `WITH_SuperLU=ON`; it links the system SuperLU library and does not alter the production solver path.

主要な build commands:

```text
git -C ../tools/elmer-hypre/src worktree add --detach ../tools/elmer-phase19-feature-gate 5a8de867068be0568f09af40fb90fee300dfbede
git -C ../tools/elmer-phase19-feature-gate apply --check docs/hypre_gpu_phase19_schur_feature.patch
cmake -S ../tools/elmer-phase19-feature-gate -B /home/symme/elmer-phase19-feature-cpu-build -G Ninja -DWITH_MPI=ON -DWITH_Mumps=ON -DWITH_Hypre=ON -DWITH_AMGX=OFF -DWITH_SuperLU=ON
cmake --build /home/symme/elmer-phase19-feature-cpu-build --parallel 4
cmake --install /home/symme/elmer-phase19-feature-cpu-build
```

前回の ParMETIS 検出失敗は、`/usr/include/parmetis` という固定仮定が、Ubuntu の `/usr/include/scotch/parmetis.h` 配置と一致しなかったことが原因だった。build helper は標準 include root から `parmetis.h` を探索する方式に変更した。新規コードにユーザー固有の WSL パスは追加していない。

## Phase B: matrix-free Schur diagnostic

対象の block sizes は `nu=84636`、`nl=2898`、`D=0`。実行ログは `results/case_p19_hypre_block_schur_diag_cpu_time5us/solver_phase19_cpu_superlu_parity.log` にある。各ベクトルについて `B^T v`、K solve、`B K^{-1}B^T v`、`Dv`、最終 `D v-B K^{-1}B^T v` の norm/min/max/nonzero/finite を記録した。全段階で `finite=T` だった。

| vector | `||B^Tv||_2` | `||K^{-1}B^Tv||_2` | `||BK^{-1}B^Tv||_2` | oracle action norm | action relative error | gate |
|---|---:|---:|---:|---:|---:|---|
| all ones | 1.71951e-08 | 8.32363e+02 | 1.06040e-07 | 1.060403652e-07 | 2.70e-07 | fail |
| alternating sign | 8.92621e-09 | 1.82733e+01 | 2.42168e-09 | 2.421675477e-09 | 5.07e-07 | fail |
| deterministic sine | 1.16776e-08 | 5.34099e+01 | 8.70324e-09 | 8.703234620e-09 | 2.56e-07 | fail |
| basis-like 17 | 1.42681e-10 | 3.06536e-01 | 3.99849e-11 | 3.998488629e-11 | 1.15e-09 | fail |

再現 command:

```text
ELMER_HOME=../tools/elmer-phase19-feature-cpu-install \
LD_LIBRARY_PATH=../tools/elmer-phase19-feature-cpu-install/lib/elmersolver:../tools/hypre-cuda-install/lib \
  ../tools/elmer-phase19-feature-cpu-install/bin/ElmerSolver_mpi generated/cases/case_p19_hypre_block_schur_diag_cpu_time5us.sif
python scripts/analysis/validate_matrix_free_schur.py --matrix case_p19_hypre_flexgmres_mgr_cpu_time5us_smoke_1step_a.dat --rows 87534 --c-start 84637 --elmer-prefix case_p19_hypre_block_schur_diag_cpu_time5us --output artifacts/hypre_phase19_schur/matrix_free_validation_rebuilt_cpu_superlu.json
```

vector 自体は 4/4 pass（最大相対誤差 `1.32e-16`）。action はノルムの絶対値と sparsity pattern が一致するが、最大相対誤差 `5.0720165e-7` で厳格値に届かない。これは UMFPACK 版でも同じ桁の差が出たため、残る問題は stage 欠落ではなく、非常に小さい Schur action に対する factorization/backend の数値差である。許容値は変更していない。

### SuperLU parity check

Elmer の diagnostic が保存した `K^-1 B^T v` と、同じ明示行列に対する SciPy `splu(..., permc_spec="COLAMD")` を比較した。K solve 相対差は all-ones `7.945e-6`、alternating `2.577e-6`、sine `3.727e-6`、basis `2.216e-12`。system SuperLU と SciPy が同一実装ではなく、同じ ordering 指定だけでは strict parity にならないことを確認した。結果は [superlu_parity_cpu.json](../artifacts/hypre_phase19_schur/superlu_parity_cpu.json) に保存した。なお、Elmer action と explicit matrix から Elmer K solve を再構成した差も `1.7e-7〜2.3e-7` であり、solver 差だけでなく block 抽出行列の fingerprint も次に固定すべきである。

### 直した根本原因

診断 loop の 1 回目後に `Aij` が B を指したままになり、2〜4 回目の `B^T v` が誤って B を使っていた。各 q の冒頭で `Aij => Submatrix(1,2)%Mat` を再設定した。また、`B^T) の行数が K の 84636 より小さい 81501 であるため、出力 RHS/work array をゼロ初期化し、未初期化 tail が混入しないようにした。

## 未達ゲートと次の一手

CPU lower/full one-step、GPU、transient は未実行・未承認。未達ゲートは「独立 SciPy/SuperLU action 相対誤差 `<=1e-10`」である。

次の一手は、block 抽出後の `K,B,B^T,D` の fingerprint を保存し、Elmer と oracle が同一行列を使うことを確認した上で、同一 SuperLU 実装・同一オプションの parity check を再実行すること。これが pass するまで production promotion は行わない。

## 変更ファイル

- [docs/hypre_gpu_phase19_schur_feature.patch](hypre_gpu_phase19_schur_feature.patch): base から clean apply できる feature patch。Schur stage diagnostics、4-vector reset、optional serial SuperLU gate を含む。
- feature source diff: `CMakeLists.txt`、`fem/src/CMakeLists.txt`、`fem/src/BlockSolve.F90`、`fem/src/SolveSuperLUStandard.c`、`fem/src/SOLVER.KEYWORDS`、`fem/src/SParIterSolver.F90`、`fem/src/SolveHypre.c`。
- [scripts/support/build_elmer_hypre_gpu_wsl.ps1](../scripts/support/build_elmer_hypre_gpu_wsl.ps1): ParMETIS header discovery を portable 化。
- [scripts/run_hypre_gpu_wsl.ps1](../scripts/run_hypre_gpu_wsl.ps1): HYPRE tag suffix、installed module/lib path、tag forwarding を修正。
- [scripts/analysis/validate_matrix_free_schur.py](../scripts/analysis/validate_matrix_free_schur.py): strict validator は変更なし。
- [artifacts/hypre_phase19_schur/matrix_free_validation_rebuilt_cpu_superlu.json](../artifacts/hypre_phase19_schur/matrix_free_validation_rebuilt_cpu_superlu.json): 再実行結果。
- [artifacts/hypre_phase19_schur/matrix_free_stage_metrics_cpu_superlu.json](../artifacts/hypre_phase19_schur/matrix_free_stage_metrics_cpu_superlu.json): 全 stage の診断要約。
- [artifacts/hypre_phase19_schur/superlu_parity_cpu.json](../artifacts/hypre_phase19_schur/superlu_parity_cpu.json): Elmer K solve と SciPy/SuperLU の parity 結果。
- `tools/elmer-phase19-feature-gate`: isolated build source worktree。通常の `main` worktree は変更していない。
