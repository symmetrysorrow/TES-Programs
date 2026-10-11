# DualTES transient performance: handoff package

## Purpose

This package is a reproducible evidence bundle for an independent performance analysis. It focuses on the user's target of bringing the DualTES transient runtime close to twice the single-pixel runtime, and on locating time that the current short-run instrumentation does not split out. No production solver setting is changed by this evidence commit.

## Key runs

| Run | Physical steps | MPI ranks | Solver elapsed (CPU, REAL) | HeatSolve Assembly sum | HeatSolve Solve sum |
|---|---:|---:|---:|---:|---:|
| DualTES 4-step diagnostic | 4 | 6 | 121.81 s | 103.80 s | 14.29 s |
| Single-pixel 4-step diagnostic | 4 | 6 | 90.26 s | 83.58 s | 2.38 s |
| DualTES pulse baseline | 225 output steps (see input SIF) | 6 | 8653.66 s (2 h 24 min 13.66 s) | 8279.85 s | 369.67 s |

The two short diagnostics use the same step sequence: 18 us once, 1 us twice, then 1 ns once. Each had six nonlinear HeatSolve iterations. Their total-runtime ratio is 1.35. This is a useful diagnostic comparison, not a precise full-run scaling estimate: the mesh families differ and the short-run restart fields have different convergence histories.

The two short runs both report Elmer 26.1, six MPI ranks, one OpenMP thread per rank, HYPRE linked, and the same HeatSolve DLL SHA256 (`be54d5d6087c38ce14bb02241e292ac9e3f2baf2f846a70543fbc39ca2b4002f`). The short-run SIFs and manifests are included. The single-pixel mesh is `mesh_hybrid_fullconf_h8`; the DualTES mesh is `mesh_toml_95333801e1a3`.

## Findings from the current evidence

### 1. The apparent assembly gap is not yet a valid subtraction

For the 4-step DualTES diagnostic, the six `HeatSolve: iter ... Assembly` values sum to 103.80 s. The Phase24 detailed event values emitted by the same log sum to:

- element traversal: 14.836 s
- material-property evaluation: 2.108 s
- local element calculation: 18.709 s
- global sparse insertion: 1.308 s
- boundary assembly: 0.243 s
- matrix-format conversion: 0 s

These detailed counters total 37.204 s, but they use `RealTime()` wall-clock intervals. The outer Assembly line is computed with `CPUTime()`. Do not label the difference (66.596 s) as measured hidden wall time.

The included `source/HeatSolve.F90` snapshot shows the Phase24 bulk interval starts after `DefaultInitialize()` and the bulk loop. The outer Assembly interval continues through `DefaultDirichletBCs()`, optional preconditioner construction, and operator bookkeeping. Those are candidate unmeasured regions; their individual costs have not yet been measured. The source snapshot is from the available `tools/elmer-hypre` tree (SHA256 `965fa34d9d051467d821e87c62806c8cea866e1478e3fc83f6bc5d20d2fe3f50`). The run manifest records the loaded HeatSolve DLL hash, but build provenance does not prove this source snapshot exactly matches that DLL. Treat it as code evidence, not a verified source-to-binary mapping.

### 2. There is a large outer-vs-inner HYPRE timing difference

In the 4-step DualTES log, the outer Phase24 linear-solve intervals sum to 14.279 s, while `SolveHypre: Solution time (method 601)` sums to 0.916 s. One outer solve interval is 12.119 s while its corresponding inner method time is 0.051 s. The single-pixel sums are 2.381 s outer versus 0.299 s inner.

In the full DualTES baseline, the 507 outer Solve intervals sum to 369.67 s, while the 507 inner `SolveHypre: Solution time` entries sum to 238.948 s. This indicates that timing only the Krylov-method line misses substantial time in at least some solve calls. The logs alone do not identify whether this is HYPRE setup, matrix transfer/conversion, synchronization, or another wrapper cost. Validate the meaning of each timer before attributing the difference.

### 3. The full transient baseline is available

The full DualTES baseline completed successfully with solver-reported CPU and REAL time both 8653.66 s. It used 225 saved time steps (507 HeatSolve nonlinear iterations). The full SIF, solver log, manifest, time-series CSVs, iteration CSVs, result plot, and final flush audit are included under `dualtes-full/`.

## Evidence map

- `dualtes-short/`: four-step DualTES profile, generated SIF, exact JSON, driver, mesh header, solver log, manifests, iteration and output series.
- `singlepixel-short/`: matching four-step single-pixel profile and the same evidence types.
- `dualtes-full/`: the previously completed 75 ms DualTES pulse run, including the plot and full log.
- `source/`: available HeatSolve source snapshot and the Phase24 instrumentation patch notes.

The manifests preserve input and runtime artifact hashes. Large mesh binaries and restart state files are deliberately omitted; their hashes and paths remain recorded in the manifests.

## Recommended next probe (still short; no full 75 ms run)

1. Rebuild an isolated diagnostic HeatSolve module from the exact source used for the run, or establish the source-to-DLL mapping first.
2. Add same-clock (`RealTime`) measurements around `DefaultInitialize()`, `DefaultDirichletBCs()`, preconditioner construction, and the remaining pre/post-bulk assembly code. Also add same-clock start/end timing around the complete outer assembly and complete outer solve.
3. Re-run only the four-step DualTES and single-pixel cases, retaining six MPI ranks and their existing settings.
4. Record per-rank and global-max timings explicitly. Avoid subtracting rank-local CPU time from wall-clock timers.
5. Determine whether HYPRE setup/reuse is already enabled and whether the 12.119 s first/second-call outlier is one-time or reproducible before recommending a default change.

## Request to the next analyst

Please independently analyze this bundle and answer:

1. Which measured interval(s) explain the DualTES vs single-pixel delta, after comparing like-for-like clocks and ranks?
2. What causes the outer HYPRE-solve time to exceed the `SolveHypre: Solution time` line, especially the 12.119 s diagnostic call and the 130.722 s aggregate gap in the full run?
3. Which minimal instrumentation or non-invasive short experiment will identify the dominant cost without running the full transient again?
4. Which optimization is most likely to move DualTES toward 2x single-pixel runtime while preserving numerical accuracy, and what evidence is needed before changing defaults?

Do not treat the 66.596 s Assembly/event difference as an established wall-time bucket. Do not infer the HYPRE setup cause solely from the log.

