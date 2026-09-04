# Phase 19 HYPRE GPU / mortar production status

更新日: 2026-09-04

## 結論

Phase A のクリーン適用と CPU HYPRE ビルドは完了した。Phase B は修正版の serial SuperLU wrapper で 4 ベクトルを再実行し、`B^T v`、K solve、`B K^{-1}B^T v`、`Dv`、最終 action を保存した。

actual Elmer block だけで作った self-consistency oracle では、stage の matvec/reconstruction は一致したが、SciPy K solve を含む strict oracle gate は 4 vector 中 3 vector で未達だった。したがって指示どおり、CPU lower/full、GPU、短時間 transient、MPI へは進めていない。raw block mismatch だけで matrix-free operator を不合格にはしていない。

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

vector 自体は 4/4 pass（最大相対誤差 `1.32e-16`）。monolithic oracle との action 差は最大 `5.0716e-7` だが、これは別の block-equivalence 問題として扱う。actual Elmer block に対しては `B^T v`、`B*ku`、`Dv`、`y=dv-bku` の emitted-stage consistency が成立した一方、SciPy K oracle を通した `B K^{-1}B^T v` は all-ones `2.51e-10`、alternating `8.93e-11`、sine `8.61e-10`、basis `3.40e-16` で、strict `1e-10` は 3 vector が未達だった。

### SuperLU parity check

actual `K_elmer` に対する Elmer emitted solve の backward residual は、相対 `3.10e-14〜5.03e-8`、componentwise `1.31e-15〜4.79e-7`。同じ `K_elmer` に対する SciPy solve の residual も併記した。system SuperLU と SciPy は同一 binary ではないため、K solution vector差は backend/pivoting 差を含む。詳細は [superlu_parity_cpu.json](../artifacts/hypre_phase19_schur/superlu_parity_cpu.json) の `self_consistency` に保存した。

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

`*_K/B/Bt/D.triplets` を logical shape（短い Bt は zero padding）へ構築し、actual blockだけで SciPy SuperLU oracleを作った。emitted stage の内部matvec検証は全4 vectorで通過したが、異なる SuperLU binary による K solve差が `B*ku` と Schur oracleへ伝播し、推奨 strict gateは未達となった。`B vs Bt^T` は相対Frobenius `0`、最大差 `0`、nnz差 `0` で PASS。Btの81501行保存は、monolithicの後続行にnonzeroがなくzero paddingが安全だった。

### 直した根本原因

診断 loop の 1 回目後に `Aij` が B を指したままになり、2〜4 回目の `B^T v` が誤って B を使っていた。各 q の冒頭で `Aij => Submatrix(1,2)%Mat` を再設定した。また、`B^T` の行数が K の 84636 より小さい 81501 であるため、出力 RHS/work array をゼロ初期化し、未初期化 tail が混入しないようにした。

## 未達ゲートと次の一手

CPU lower/full one-step、GPU、transient は未実行・未承認。未達ゲートは actual block oracle の strict Schur self error `<=1e-10` と、backend差を含む K solve評価である。

今回の次段階として actual-block self-consistency、`B vs Bt^T`、raw/canonical block equivalence、monolithic perturbation、scale-aware diagonal threshold を保存した。次は K backend差を含まない同一direct-solver oracle、または同一 binaryでの比較を用いて strict self gateを再確認すること。これが passするまで production promotionは行わない。

## 変更ファイル

- [docs/hypre_gpu_phase19_schur_feature.patch](hypre_gpu_phase19_schur_feature.patch): base から clean apply できる feature patch。Schur stage diagnostics、4-vector reset、optional serial SuperLU gate を含む。
- feature source diff: `CMakeLists.txt`、`fem/src/CMakeLists.txt`、`fem/src/BlockSolve.F90`、`fem/src/SolveSuperLUStandard.c`、`fem/src/SOLVER.KEYWORDS`、`fem/src/SParIterSolver.F90`、`fem/src/SolveHypre.c`。
- [scripts/support/build_elmer_hypre_gpu_wsl.ps1](../scripts/support/build_elmer_hypre_gpu_wsl.ps1): ParMETIS header discovery を portable 化。
- [scripts/run_hypre_gpu_wsl.ps1](../scripts/run_hypre_gpu_wsl.ps1): HYPRE tag suffix、installed module/lib path、tag forwarding を修正。
- [scripts/analysis/validate_matrix_free_schur.py](../scripts/analysis/validate_matrix_free_schur.py): strict validator は変更なし。
- [artifacts/hypre_phase19_schur/matrix_free_validation_current.json](../artifacts/hypre_phase19_schur/matrix_free_validation_current.json): 修正版 wrapper 後の 4-vector 再実行結果。
- [artifacts/hypre_phase19_schur/matrix_free_stage_metrics_cpu_superlu.json](../artifacts/hypre_phase19_schur/matrix_free_stage_metrics_cpu_superlu.json): 全 stage の診断要約。
- [artifacts/hypre_phase19_schur/superlu_parity_cpu.json](../artifacts/hypre_phase19_schur/superlu_parity_cpu.json): Elmer K solve と SciPy/SuperLU の parity 結果。
- raw `*_K/B/Bt/D.triplets` と `_vN/_btN/_kuN/_bkuN/_dvN/_yN.dat` は validation中の外部artifactとして `results/case_p19_hypre_block_schur_diag_cpu_time5us/` に保持している。長期Git管理には含めない。
- `tools/elmer-phase19-feature-gate`: isolated build source worktree。通常の `main` worktree は変更していない。

## Explicit verdict

- SuperLU wrapper: FIXED（CSR→CSC、caller array 非破壊、合成非対称テスト PASS）
- Matrix-free implementation self-consistency: FAIL（emitted matvec/reconstructionはPASS、SciPy actual-block oracleは3/4 strict FAIL）
- B vs Bt^T: PASS（相対Frobenius 0、最大差0、nnz差0）
- Block extraction vs monolithic: NUMERICALLY CLOSE（raw SHA/shapeは不一致、common entry exact、extraはmachine-epsilon級）
- Matrix-free Schur implementation: INVALID（raw SHA mismatchではなく、actual-block oracle strict gate未達）
- Lower CPU: NOT RUN
- Full CPU: NOT RUN
- GPU ready: NO
