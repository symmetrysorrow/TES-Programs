# Phase24 short difference campaign

Generated: `2026-09-15T15:32:05.363232+00:00`
Project: `D:\Github\TES-Programs\Elmer-Projects\artifacts\phase24_short_difference_campaign\phase24_short_campaign.json`

## Decision rule

- If BDF2 reuse-OFF moves materially toward BDF1 reuse-OFF, reuse/lifecycle is implicated.
- If BDF1 reuse-OFF matches the CPU reference while BDF2 reuse-OFF does not, BDF2/history handling is implicated.
- If both HYPRE variants differ similarly from CPU/MUMPS, investigate Phase24 assembly, HYPRE tolerance/conditioning, or the runtime UDF before changing the time integrator.

The CPU/MUMPS reference is not a bit-identical solver stack: it uses the historical Phase23 binary and BDF1. The campaign therefore localizes the cause; it is not itself a final numerical qualification.

## Runs

| variant | exit code | comparison |
|---|---:|---|
| BDF2 / HYPRE reuse ON | 0 | `D:\Github\TES-Programs\Elmer-Projects\artifacts\comparison\phase24_bdf2_reuse_on_1us` |
| BDF2 / reuse OFF | 0 | `D:\Github\TES-Programs\Elmer-Projects\artifacts\comparison\phase24_bdf2_reuse_off_1us` |
| BDF1 / reuse OFF | 0 | `D:\Github\TES-Programs\Elmer-Projects\artifacts\comparison\phase24_bdf1_reuse_off_1us` |

All solver output is retained in the corresponding `results/case_phase24_short_*_1us/` directory. Each comparison directory contains `summary.md`, `summary.json`, `aligned_difference.csv`, `checkpoints.csv`, and `comparison.png` when plotting is available.
