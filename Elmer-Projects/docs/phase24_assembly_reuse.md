# Phase24 assembly/reuse implementation baseline

Phase24 の production KPI は solver iteration 数ではなく、同一 TES physics・同一
case の transient end-to-end wall time とする。GPU backend は保持するが、最初の
production candidate は CPU conformal HYPRE とする。

## 現状監査

現行の対応する Elmer source tree は `tools/elmer-hypre` にあり、このリポジトリの
Git管理対象外である。監査時点では次の構造が確認できる。

- `HeatSolve.F90` の legacy solver は bulk element loop がシリアルで、各 element の
  最後に `DefaultUpdateEquations` を呼ぶ。
- legacy call は `CRS_GlueLocalMatrix` に入り、global row/column から CRS row 内を
  毎回探索する。
- Elmer には `HeatSolveVec.F90` と `CRS_GlueLocalMatrixVec` が既にあり、vectorized
  glue、element coloring、OpenMP loop の実装がある。ただしこれは TES の legacy
  circuit/特殊BCを含む production route と同値とは限らないため、そのまま backend
  切替には使わない。
- `HeatSolve.F90` は現在も `CPUTime` を assembly/solve の境界に使っている箇所が
  ある。CPU time と wall time を production acceptance に混在させない。

したがって最初の変更は、physics を変えずに計測可能な optional path とする。
legacy source への適用点は次の順序で固定する。

1. `HeatSolve.F90` に optional な Phase24 wall events を追加する。粗い確認用に
   UDF/circuit、bulk assembly、boundary assembly、matrix finalization、linear solve の
   区間を出し、artifact用にはbulkを material-property evaluation、local element
   calculation、global sparse insertion、残余のelement traversalへ分けてREAL timeで
   出力する。親bulk区間はartifactへ二重計上しない。
2. `DefaultUpdateEquations(..., VecAssembly=.TRUE.)` を明示的な opt-in として legacy
   route で比較する。これは insertion destination reuse の候補を測るための段階で、
   matrix value reuse ではない。
3. correctness gate を通過した場合のみ、`CRS_GlueLocalMatrixVec` に element-index
   keyed insertion-position cache を追加する。cache key は matrix structure hash、
   element id、local dof layout、permutation signature を含める。
4. geometry/shape-gradient cache は P1 tetra のみに限定し、coordinate-dependent
   UDF、moving mesh、radiation、phase-change derivative がある場合は無効化する。
5. CPU element parallelism は、まず local calculation の thread-local buffers と
   serial merge で測る。global insertion の並列化は existing coloring の validation
   後に行い、大規模 critical section は使わない。

## 実装済みの第一段階

`tools/elmer-hypre/src/fem/src/CRSMatrix.F90` の Phase24 cache は、scalar P1 tetra の
4 x 4 local matrixについて、初回だけCSR destinationを解決し、以後は16個の値を
直接加算する。cacheは element connectivity、matrix object、CSR `Rows`/`Cols` array
target、配列サイズをguardに持つ。CSR配列が再確保された場合は古いposition mapを
再利用しない。

`HeatSolve.F90` のopt-in vector assemblyは、3D scalar P1 tetra meshであることを
初期化時に確認する。対象外のmeshではPhase24 pathを自動的に無効化し、generic
Elmer assemblyへ戻る。デフォルトopt-in flagは引き続きfalseである。

同じopt-in pathではbody単位のconstant scalar `Heat Capacity` と `Heat Conductivity`
を一度だけ解決する。MATC/UDF、異方性、その他のdynamic propertyはcache対象外で、
従来のgeneric lookupを使う。

## Artifact contract

`phase24_profile.py` が生成する wall profile は、nested event の子区間を親から
引いて additive にする。これにより、次が成立する。

`sum(measured wall components) + unclassified ~= total solver wall`

CPU-only event は拒否する。以下のcorrectness gateを満たさないcandidateは、速度に
関係なくrejectする。

- matrix sparsity pattern
- matrix values
- RHS
- full temperature field
- TES volume-average / absorber temperature
- current / resistance / Joule power
- pulse response
- solver convergence

実行例:

```powershell
python scripts/analysis/phase24_profile.py profile `
  --events phase24_events.json --total-wall 12.5 --case production_smoke `
  --output artifacts/assembly_wall_profile.json

python scripts/analysis/phase24_profile.py log `
  --log results/production_smoke/solver.log --total-wall 12.5 --case production_smoke `
  --output artifacts/assembly_wall_profile.json

python scripts/analysis/phase24_profile.py validate `
  --baseline baseline_gate.json --candidate cached_gate.json `
  --output artifacts/phase24_regression_gate.json

python scripts/analysis/phase24_profile.py classify `
  --baseline-profile artifacts/assembly_wall_profile.json `
  --candidate-profile artifacts/cached_profile.json `
  --gate artifacts/phase24_regression_gate.json `
  --output artifacts/phase24_backend_recommendation.json
```

実測前に benchmark artifact を作成して「改善した」と扱わない。`unclassified` が
大きい profile は instrumentation incomplete として扱い、まず Step A を完了する。

## 実装順と停止条件

Step A の wall breakdown → insertion destination cache → geometry metadata → constant
material metadata → local OpenMP → colored global assembly → static/dynamic split →
HYPRE structure allocation reuse の順とする。各段階で total wall、assembly wall、
insertion、HYPRE setup/solve、peak memory、全 correctness gate を保存する。

同一caseで改善が 2% 以下なら `NO_MEASURABLE_GAIN` 相当として次段階へ進み、複雑化
だけを残さない。`Linear System Refactorize = False` のような blind timestep reuse は
この段階では禁止する。

## Build unblock and Stage 2 evidence

The existing configured native build uses MSYS2 UCRT64 GNU compiler/MPI wrappers with Microsoft MPI import libraries. Its CMake cache already points at `C:/msys64/ucrt64/bin/mpif90.exe`, `mpicc.exe`, `mpicxx.exe`, `libmsmpi.dll.a`, and `C:/Program Files/Microsoft MPI/Bin/mpiexec.exe`. The apparent MPI/compiler failure was caused by PowerShell missing `C:/msys64/ucrt64/bin` and `C:/msys64/usr/bin` on `PATH`; the Fortran frontend exited with Windows status `0xC0000139`. Prepending those directories restored compilation using the existing cache. No MPI implementation was installed or switched, and WSL OpenMPI was not mixed into the native build.

Stage 2 now retains fixed P1 tetra geometry and compact static/dynamic/generic lists in `HeatSolve.F90`, computes reusable static stiffness/spatial-mass buffers, and performs static local assembly with OpenMP plus cached CSR insertion and atomics. The generic path remains authoritative for the production hybrid prism/tetra mesh. The bounded CPU HYPRE smoke improved solver wall time from 20.75 s to 11.18 s and assembly from 13.51 s to 2.23 s. The smoke and gate are recorded in `phase24_production_before_after.json` and `phase24_correctness_gate.json`; the design and fallback rules are in `phase24_fast_path_design.md`.
