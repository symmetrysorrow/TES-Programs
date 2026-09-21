# Phase24 refined-reference causal validation

実行日: 2026-09-22 (Asia/Tokyo)

## 判定

`--refine-reference` の MUMPS steady refinement と、同じ refined checkpoint
からの MUMPS/HYPRE transient は Full ElmerSolver で正常完了した。
refined reference の最初のBDF1 outer solveは、従来の `-0.76865 mK` jumpではなく
`+0.113118 µK` (`nl1` before → outer-after) だった。

従って first-step artifact の初期化依存性は実験的に再現され、jumpは少なくとも
桁違いに縮小した。ただし refinement は steady fixed point 自体を大きく移動させた
ため、COMSOL parityを保ったまま「Stycast primal residualだけを除去した」とは
まだ分離できない。post-refinement の full matrix residual capture は今回の
production runでは有効化していない。

## Exact command

```powershell
python scripts/support/run_phase24_gate4_5_nomortar.py `
  --window 40us `
  --backend both `
  --mortar `
  --refine-reference `
  --base-project artifacts/phase24_gate3_hypre_steady_s32m2/phase24_gate3_hypre_steady.json `
  --mesh mesh_singlepixel_gpu_fine_stycast32_mortar `
  --reference-case case_phase24_g3_s32m_s32m2_10f17b
```

実行時のsolver commandはStage11 `ElmerSolver.exe`、MPI procs=1、
`--runtime-bin D:\Github\TES-Programs\tools\elmer-hypre\install-stage11\bin`、
`--toolchain-bin C:\msys64\ucrt64\bin`。`ELMER_HOME` は同Stage11 prefix、
`ELMER_LIB` は未設定だった。child PATH の先頭は
`C:\msys64\ucrt64\bin;D:\Github\TES-Programs\tools\elmer-hypre\install-stage11\bin;D:\Github\TES-Programs\tools\elmer-hypre\install-stage11\share\elmersolver\lib`。

runtime hashes:

- ElmerSolver: `64922c6d98390e7edebb16db87df8118b46fda4b68ed9a483eb7b8e556a26448`
- libelmersolver: `b50a6a98a824bf5c16d9e33ce8b504b3d1b5ca9f286250722b28dfc465fecfaa`
- HeatSolve: `851efd3f28d823d6d6e72958509233bec1382c1376787b76d537e0092df0968f`

## Steady refinement

| quantity | Gate3 HYPRE restart | refined MUMPS fixed point |
|---|---:|---:|
| TES temperature | 168.569115 mK (checkpoint; Gate3 summary 168.569130 mK) | 166.556911 mK |
| current | 143.568115 µA (checkpoint; Gate3 summary 143.567589 µA) | 218.667666 µA |
| resistance | 15.522836 mΩ | 8.852228 mΩ |
| power | 319.906331 pW | 423.274148 pW raw / 423.228618 pW relaxed |

Refinement converged in 84 nonlinear iterations; final `ComputeChange` relative
change was `9.893626e-9`; solver log contains `MAIN: *** Elmer Solver: ALL DONE ***`.
The shifts are `-2.012204 mK` and `+75.099551 µA`, so the original Gate3 COMSOL
current parity (`+0.358282%` in the Gate3 summary) does not carry over.

## `.result` / `.state` handoff

The refined state contains, in order, `T`, `Current`, `Resistance`, `Power`,
`PreviousCurrent`:

```text
0.16655691141885667  0.00021866766560428157  0.008852228329204792
4.232286181467137e-10  0.00021866766560428157
```

These values match the final refinement iteration's T/current/resistance and
relaxed power; `PreviousCurrent` is committed to the refined current. The refined
`.result` and `.state` are separate from the Gate3 inputs.

## First-step diagnostic

The production 40-us run includes the requested pulse-OFF BDF1 `dt=18 µs` first
step. The first timestep values were identical for both backends:

- `nl1 x_before`: `166.556911419 mK`
- `nl1 outer_after`: `166.557024537 mK`
- `ΔT`: `+0.113118 µK`
- first accepted pre-pulse output: `166.557024537 mK`, `218.663253 µA` (HYPRE)
- first accepted post-pulse output at `+0.001 µs`: MUMPS `166.557079389 mK / 218.662000 µA`; HYPRE `166.557024558 mK / 218.663230 µA`
- first-step circuit residual: `4.141108e-14 W` at nl1→outer-after

Compared with `-768.65 µK`, the first-step temperature change is smaller by
approximately 6.8 million times in magnitude. Numeric post-refinement primal,
constraint, worst-interface-row, and componentwise backward-error values are
`not recorded` because this production run did not emit the full restriction
matrix capture. The refinement nonlinear/circuit residual above must not be
mislabelled as the assembled primal residual.

## 40-us transient completion and parity

Both transient solvers completed normally (`exit_code=0`, `ALL DONE`). The
generated time grid ended at `38.751 µs` after the pulse, not 40.000 µs; the
runner therefore exited while its 40-us comparison guard was evaluating. The
offline comparison below uses the common available endpoint and records this
coverage limitation explicitly.

| comparison | max difference | RMSE | result |
|---|---:|---:|---|
| MUMPS vs COMSOL, baseline-subtracted | 0.129232 µA | 0.054824 µA | target fail |
| HYPRE vs COMSOL, baseline-subtracted | 0.616476 µA | 0.216891 µA | target fail |
| HYPRE vs MUMPS | 0.001167 µA | 0.001082 µA | backend close |

Absolute baseline parity is not acceptable after refinement:

- COMSOL: `143.055049 µA`
- MUMPS baseline: `218.663013 µA` (`+75.607964 µA`, `+52.852356%`)
- HYPRE baseline: `218.663242 µA` (`+75.608193 µA`, `+52.852516%`)

Thus the baseline-subtracted waveform and absolute baseline must not be conflated.

## Input protection and provenance

The original Gate3 `.result` and `.state` SHA256 values were unchanged:

- result before/after: `3bb73ac319f7302eb3a9fb6dee2a8c558a7a7a54426671d874c7e31358948424`
- state before/after: `9c2eca9bbb26fc6fab1949511787b87350ccdc370658ed2a41d458959f5971a6`

Relevant logs/manifests are under `results/`:

- `case_phase24_restart_refine_fine_stycast32_mortar_mumps_mortar`
- `case_phase24_g45_fine_stycast32_mortar_40us_mumps_mortar`
- `case_phase24_g45_fine_stycast32_mortar_40us_hypre_mortar`

The runner's generated project, refinement SIF, transient SIFs, manifests, logs,
iteration CSVs, and waveform CSVs retain the exact input/output paths and hashes.

## Final causal interpretation

Case A's observable criterion passes: refining the restart removes the known
first-step jump. The stronger claim that the entire `~0.769 mK` artifact was
uniquely the Gate3 Stycast interface primal residual remains provisional because
the refined steady state moved by 2.012 mK and 75.10 µA, and no post-refinement
full saddle residual capture was collected. The next task is a read-only full
restriction capture of the refined steady solve/restart followed by a corrected
40-us endpoint comparison; physics, production SIF, solver formulation,
preconditioner, scaling, and tolerances remain unchanged.
