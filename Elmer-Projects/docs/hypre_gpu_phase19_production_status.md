# Phase 19 HYPRE GPU / mortar production status

更新日: 2026-09-04

## 結論

Phase A のクリーン適用と CPU HYPRE ビルドは完了した。Phase B は修正版の serial SuperLU wrapper で 4 ベクトルを再実行し、`B^T v`、K solve、`B K^{-1}B^T v`、`Dv`、最終 action を保存した。

actual Elmer block だけで作った self-consistency oracle では、stage の matvec/reconstruction は一致したが、SciPy K solve を含む strict oracle gate は 4 vector 中 3 vector で未達だった。これは cross-backend sensitivity として保存し、implementation correctness の gate から分離した。

Phase C では、production `BlockSchurActionWork` の short-Bt tail zeroing を修正し、同じ `elmersolver` にリンクした `BlockSchurSuperLUOracle` で4 RHSを再検証した。same-binary K parity と same-binary Schur parity は4/4 PASS、stage algebraもPASSだった。Schur diagonalは `scale=4.59957e-11`、`threshold=6.85389e-19`、`min abs=4.68324e-13`、`max abs=4.59957e-11`、`min/max ratio=1.01819e-2`、negative `2898`、positive `0`、near-zero `0` である。

SciPy比較は診断用に残す。`Schur <= 1e-10` は2/4 fail、`B*ku <= 1e-12`を含むcomposite cross-backend gateは3/4 failであり、これをsame-binary implementation gateとは混同しない。lower は同一条件の CPU baseline を実行したが約45分経過時点で outer solve が終了せず、`INCOMPLETE` として不採用にした。full は lower が受理可能な結果を出していないため未実行であり、GPU・MPI・transient・tuningも未実施である。

## Phase A: clean base / build

- source base: `5a8de867068be0568f09af40fb90fee300dfbede`
- feature patch: [hypre_gpu_phase19_schur_feature.patch](hypre_gpu_phase19_schur_feature.patch)
- `git apply --check`: clean
- build: isolated source worktree + Ninja, MPI/HYPRE/MUMPS/UMFPACK enabled
- HYPRE library: `tools/hypre-cuda-install/lib/libHYPRE.so` (CPU execution mode)
- executable: `tools/elmer-phase19-feature-cpu-install/bin/ElmerSolver_mpi`
- optional serial SuperLU diagnostic wrapper is enabled by `WITH_SuperLU=ON`; it links the system SuperLU library and does not alter the production solver path.
- CSR-to-CSC conversion is explicit, range-checked, and uses temporary arrays; Elmer `Rows`, `Cols`, and `Values` are not mutated.
- `BlockSchurSuperLUStandardTest` is a non-symmetric synthetic regression test and passes (`K` residual `0`, `K^T` residual `10.198039...`).
- Schur diagonal near-zero handling is scale-aware: `max(100*TINY, sqrt(epsilon)*max|DiagS|)`, with sign and threshold diagnostics logged.

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

対象の block sizes は `nu=84636`、`nl=2898`、`D=0`。実行ログは `results/case_p19_hypre_block_schur_diag_cpu_time5us/solver.log` にある。各ベクトルについて `B^T v`、K solve、`B K^{-1}B^T v`、`Dv`、最終 `D v-B K^{-1}B^T v` の norm/min/max/nonzero/finite を記録した。全段階で `finite=T` だった。

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
python scripts/analysis/validate_matrix_free_schur.py --matrix case_p19_hypre_flexgmres_mgr_cpu_time5us_smoke_1step_a.dat --rows 87534 --c-start 84637 --elmer-prefix results/case_p19_hypre_block_schur_diag_cpu_time5us/case_p19_hypre_block_schur_diag_cpu_time5us --output artifacts/hypre_phase19_schur/matrix_free_validation_current.json
python scripts/analysis/check_superlu_parity.py --matrix case_p19_hypre_flexgmres_mgr_cpu_time5us_smoke_1step_a.dat --rows 87534 --c-start 84637 --elmer-prefix results/case_p19_hypre_block_schur_diag_cpu_time5us/case_p19_hypre_block_schur_diag_cpu_time5us --output artifacts/hypre_phase19_schur/superlu_parity_cpu.json
```

vector 自体は 4/4 pass（最大相対誤差 `1.32e-16`）。以下の action error は monolithic/SciPy cross-backend 診断値であり、same-binary gateとは別である。monolithic oracle との action 差は最大 `5.0716e-7`。actual Elmer block に対して SciPy K oracle を通した `B K^{-1}B^T v` は all-ones `2.51e-10`、alternating `8.93e-11`、sine `8.61e-10`、basis `3.40e-16` で、strict `1e-10` は 2/4 pass（2/4 fail）だった。

### SuperLU parity check

actual `K_elmer` に対する Elmer emitted solve の backward residual は、相対 `3.10e-14〜5.03e-8`、componentwise `1.31e-15〜4.79e-7`。同じ `K_elmer` に対する SciPy solve の residual も併記した。system SuperLU と SciPy は同一 binary ではないため、K solution vector差は backend/pivoting 差を含む。詳細は [same_binary_parity.json](../artifacts/hypre_phase19_schur/same_binary_parity.json) の same-binary/cross-backend 分離結果に保存した。

same-binary oracle は、Elmer diagnostic と同じ `block_schur_superlu_solve`（同じ SuperLU binary、CSR→CSC wrapper、`COLAMD`）を使用した。4 vectorの solution relative/absolute error は保存精度で `0/0`、`B*ku` relative error は `1.79e-16〜3.45e-16`、Schur relative error は `1.79e-16〜3.45e-16` だった。K backward residual relative は `3.10e-14〜5.03e-8`、componentwise backward error は `1.31e-15〜4.79e-7`。結果は [same_binary_parity.json](../artifacts/hypre_phase19_schur/same_binary_parity.json) に保存した。

same-binary K parity:

| vector | solution rel. error | solution abs. error | K backward rel. | componentwise backward | gate |
|---|---:|---:|---:|---:|---|
| all ones | 0 | 0 | 5.03e-8 | 4.79e-7 | PASS |
| alternating | 0 | 0 | 8.09e-10 | 4.51e-7 | PASS |
| sine | 0 | 0 | 2.10e-8 | 5.03e-8 | PASS |
| basis | 0 | 0 | 3.10e-14 | 1.31e-15 | PASS |

same-binary Schur parity (`B*ku` and final `y`):

| vector | `B*ku` rel. error | Schur rel. error | gate |
|---|---:|---:|---|
| all ones | 1.79e-16 | 1.79e-16 | PASS |
| alternating | 1.87e-16 | 1.87e-16 | PASS |
| sine | 1.91e-16 | 1.91e-16 | PASS |
| basis | 3.45e-16 | 3.45e-16 | PASS |

same-binary oracle の再現手順（oracle生成はWSL、集計はWindows Python）:

```text
cmake --build /home/symme/elmer-phase19-feature-cpu-build --target BlockSchurSuperLUOracle --parallel 4
wsl.exe -d Ubuntu -- bash -lc 'export LD_LIBRARY_PATH=/home/symme/elmer-phase19-feature-cpu-build/fem/src:/usr/lib/x86_64-linux-gnu; .../BlockSchurSuperLUOracle <K.triplets> <btN.dat> <same_binary_kuN.dat>'
python scripts/analysis/check_superlu_parity.py --matrix case_p19_hypre_flexgmres_mgr_cpu_time5us_smoke_1step_a.dat --rows 87534 --c-start 84637 --elmer-prefix results/case_p19_hypre_block_schur_diag_cpu_time5us/case_p19_hypre_block_schur_diag_cpu_time5us --same-binary-solution-prefix results/case_p19_hypre_block_schur_diag_cpu_time5us/case_p19_hypre_block_schur_diag_cpu_time5us_same_binary --output artifacts/hypre_phase19_schur/same_binary_parity.json
```

### Block fingerprint / exact diff

Elmer が diagnostic で抽出した block と、保存済み monolithic matrix の block を canonical triplet と SHA-256 で比較した。raw fingerprint の厳密一致は `D` のみで、`K`、`B`、`B^T` は MISMATCH である。

| block | monolithic shape / nnz | Elmer shape / nnz | raw SHA match | 差分の要点 |
|---|---:|---:|---|---|
| K | 84636×84636 / 1340068 | 84636×84636 / 1449892 | no | explicit zero を含み、追加 nonzero 746、最大 `1.17e-16` |
| B | 2898×84636 / 53564 | 2898×84636 / 54440 | no | 追加 nonzero 876、最大 `2.22e-16` |
| B^T | 84636×2898 / 53564 | 81501×2898 / 54440 | no | stored row shape が短く、追加 nonzero 876、最大 `2.22e-16` |
| D | 2898×2898 / 0 | 2898×2898 / 0 | yes | exact empty block |

raw SHA/shapeは一致しないが、K/B/Bᵀの共通entry差は0で、追加entryは機械精度級。relative Frobenius difference は K `2.39e-16`、B/Bᵀ `3.02e-7`、monolithicとのSchur action差は絶対 `2.87e-14〜4.58e-20`、相対最大 `5.07e-7` だった。従って分類は `NUMERICALLY CLOSE` とし、raw不一致と作用差を別々に記録した。

### Actual Elmer block self-consistency

`*_K/B/Bt/D.triplets` を logical shape（短い Bt は zero padding）へ構築した。emitted stage の内部matvec検証は全4 vectorで通過し、same-binary oracleでも K solve と Schur action が全4 vectorで一致した。SciPy comparison は異なる SuperLU binary による backend sensitivity として分離している。`B vs Bt^T` は相対Frobenius `0`、最大差 `0`、nnz差 `0` で PASS。Btの81501行保存は、monolithicの後続行にnonzeroがなくzero paddingが安全だった。

### 直した根本原因

診断 loop の 1 回目後に `Aij` が B を指したままになり、2〜4 回目の `B^T v` が誤って B を使っていた。各 q の冒頭で `Aij => Submatrix(1,2)%Mat` を再設定した。また、`B^T` の行数が K の 84636 より小さい 81501 であるため、出力 RHS/work array をゼロ初期化し、未初期化 tail が混入しないようにした。

## 未達ゲートと次の一手

CPU lower/full one-step、GPU、transient は未承認。same-binary correctness gate は通過済みだが、lower baseline は `INCOMPLETE`、full は未実行である。未達なのは cross-backend sensitivity と production convergence/performance readiness であり、同一 binary の matrix-free correctness gate ではない。

今回の次段階として actual-block stage algebra、same-binary K/Schur parity、`B vs Bt^T`、raw/canonical block equivalence、monolithic perturbation、scale-aware diagonal threshold を保存した。production promotion は lower/full の finite な convergence baseline が得られるまで行わない。

## 変更ファイル

- [docs/hypre_gpu_phase19_schur_feature.patch](hypre_gpu_phase19_schur_feature.patch): base から clean apply できる feature patch。Schur stage diagnostics、production work-array zeroing、short-Bt回帰、optional serial SuperLU gate/oracle を含む。
- feature source diff: `CMakeLists.txt`、`fem/src/CMakeLists.txt`、`fem/src/BlockSolve.F90`、`fem/src/SolveSuperLUStandard.c`、`fem/src/SOLVER.KEYWORDS`、`fem/src/SParIterSolver.F90`、`fem/src/SolveHypre.c`。
- [scripts/support/build_elmer_hypre_gpu_wsl.ps1](../scripts/support/build_elmer_hypre_gpu_wsl.ps1): ParMETIS header discovery を portable 化。
- [scripts/run_hypre_gpu_wsl.ps1](../scripts/run_hypre_gpu_wsl.ps1): HYPRE tag suffix、installed module/lib path、tag forwarding を修正。
- [scripts/analysis/validate_matrix_free_schur.py](../scripts/analysis/validate_matrix_free_schur.py): SciPy cross-backend validator。
- [scripts/analysis/check_superlu_parity.py](../scripts/analysis/check_superlu_parity.py): stage/same-binary/cross-backend/block/monolithicを分離したparity集計。
- [scripts/analysis/short_bt_regression.py](../scripts/analysis/short_bt_regression.py): finite stale-tail regression。
- [artifacts/hypre_phase19_schur/matrix_free_validation_current.json](../artifacts/hypre_phase19_schur/matrix_free_validation_current.json): 修正版 wrapper 後の 4-vector 再実行結果。
- [artifacts/hypre_phase19_schur/matrix_free_stage_metrics_cpu_superlu.json](../artifacts/hypre_phase19_schur/matrix_free_stage_metrics_cpu_superlu.json): 全 stage の診断要約。
- [artifacts/hypre_phase19_schur/same_binary_parity.json](../artifacts/hypre_phase19_schur/same_binary_parity.json): same-binary K/Schur oracleとcross-backend比較。
- [artifacts/hypre_phase19_schur/superlu_cross_backend.json](../artifacts/hypre_phase19_schur/superlu_cross_backend.json): SciPy comparisonをdiagnostic-onlyとして分離。
- [artifacts/hypre_phase19_schur/short_bt_regression.json](../artifacts/hypre_phase19_schur/short_bt_regression.json): finite stale-tail regression結果。
- [artifacts/hypre_phase19_schur/full_cpu_one_step.json](../artifacts/hypre_phase19_schur/full_cpu_one_step.json): lower baseline 未受理のため full CPU を未実行とした記録。
- raw `*_K/B/Bt/D.triplets` と `_vN/_btN/_kuN/_bkuN/_dvN/_yN.dat` は validation中の外部artifactとして `results/case_p19_hypre_block_schur_diag_cpu_time5us/` に保持している。長期Git管理には含めない。
- `tools/elmer-phase19-feature-gate`: isolated build source worktree。通常の `main` worktree は変更していない。

## Explicit verdict

- CSR->CSC wrapper: PASS（caller array非破壊、合成非対称テストPASS）
- Short-Bt production handling: PASS（finite stale tail `123/-456`をzero化）
- B vs Bt^T: PASS（相対Frobenius 0、最大差0、nnz差0）
- Matrix-free emitted-stage algebra: PASS（4/4）
- Same-binary exact-K parity: PASS（4/4）
- Same-binary Schur parity: PASS（4/4）
- Cross-backend SuperLU sensitivity: ACCEPTABLE（diagnostic only; Schur gate 2/4 fail、composite 3/4 fail）
- Block extraction vs monolithic: NUMERICALLY CLOSE（raw SHA/shapeは不一致、common entry exact、extraはmachine-epsilon級）
- Matrix-free Schur implementation: VALID
- Lower CPU: NOT RUN（実行試行は `INCOMPLETE`、partial iterate は不採用）
- Full CPU: NOT RUN
- Ready for GPU: NO

## Phase20 handoff

The Phase19 correctness conclusion is retained, while convergence and
performance are now investigated with bounded peer probes.  Lower CPU's
45-minute incomplete run is evidence of nested-Krylov cost and stagnation,
not evidence against GPU plumbing.  Full CPU/GPU are generated and compared
under the same conditions rather than gated on lower CPU acceptance.  The
Phase20 probe and acceptance policy are documented in
[hypre_gpu_phase20_status.md](hypre_gpu_phase20_status.md).
