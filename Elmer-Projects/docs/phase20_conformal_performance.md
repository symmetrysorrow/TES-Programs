# Phase20 conformal/shared-node HYPRE GPU performance

This report continues commit `b91d6bcd` after the correctness gates were
closed.  It measures the validated conformal route with the same production
physics, `1e-8` HYPRE tolerance, and one MPI rank.  It does not promote the
GPU route to a new correctness claim; the existing Phase20 parity gates remain
the acceptance basis.

## Decision

`GPU_SOLVER_EFFECTIVE_BUT_HOST_BOUND`

The GPU backend improves the linear-system measurement as the mesh grows, but
the production transient does not yet show a clear end-to-end benefit.  The
next target is host-side assembly/integration and a correctness-guarded matrix
update path, not an unconstrained AMG parameter sweep.

## Steady size crossover

| mesh | nodes / tetrahedra | GPU linear speedup | GPU wall speedup |
|---|---:|---:|---:|
| small | 26,705 / 132,458 | 0.40x | 0.79x |
| production | 90,872 / 459,683 | 1.28x median | 1.03x median |
| medium | 225,821 / 1,170,252 | 2.90x | 1.15x |
| fine | 382,507 / 1,978,314 | 4.66x | 1.28x |

The 90k point used three CPU and three GPU runs.  Median wall times were
17.38 s CPU and 16.93 s GPU; the observed run-to-run standard deviations were
0.38 s and 0.08 s.  CPU/GPU field parity passed at all four sizes.

## Transient prefixes

| prefix | CPU wall | GPU wall | GPU speedup |
|---|---:|---:|---:|
| 7 steps | 38.93 s | 38.86 s | 1.002x |
| 50 steps | 205.50 s | 213.59 s | 0.962x |

The 7-step and 50-step CPU/GPU fields and electrical series passed the
transient parity checks.  The 50-step result is especially informative:
solver-only acceleration does not translate into end-to-end acceleration on
this current production prefix.

## Breakdown and setup/reuse

The runner now records `WALL_SECONDS` around the complete solver invocation and
output collection.  Logs retain both Elmer's CPU/REAL totals and the native
HYPRE setup and solution timers.  The observed decomposition is:

`wall = outside-solver wall + HYPRE setup + HYPRE solve + host solver remainder`

The host remainder is deliberately labeled as a combined bucket: current
runtime logs do not independently timestamp FEM assembly, circuit/UDF,
matrix conversion, and result I/O.  At 90k, the host remainder is about 11.8–
12.2 s per steady run while HYPRE setup plus solve is about 4.5 s.  In the
50-step prefix, HYPRE setup totals 19.10 s CPU versus 26.57 s GPU, while total
wall is 205.50 s versus 213.59 s.

The transient logs show HYPRE setup on every nonlinear solve.  Explicit
`Linear System Refactorize = False` probes were rejected: both CPU and GPU
probes failed to converge at timestep 2.  The temperature-dependent matrix
cannot use unconditional preconditioner reuse.  A future `SolveHypre3`-style
matrix update must prove convergence and physical parity before reuse is
enabled.

GPU monitoring during a representative 90k run measured mean SM utilization
11.5% and maximum 74%, with mean memory-controller utilization 22.1%.  This
burst pattern is consistent with a host-bound workflow whose GPU work is
concentrated in HYPRE setup/solve windows.  Nsight Systems capture was
attempted but did not finalize under the current WSL runtime; no profiler file
was retained.

## CPU Mortar replacement baseline

The retained validated CPU Mortar reference reports 35.55 s solver REAL time,
but it uses a different mesh/runtime and is not an apples-to-apples production
transient.  Therefore no strict CPU Mortar → conformal GPU replacement
speedup is claimed.  A common-mesh Mortar transient is the required follow-up
for that number.

Machine-readable evidence:

- `artifacts/phase20_conformal/phase20_performance_acceptance.json`
- `artifacts/phase20_conformal/wall_time_breakdown.json`
- `artifacts/phase20_conformal/phase20_performance_crossover.json`
- `artifacts/phase20_conformal/phase20_transient_performance.json`
- `artifacts/phase20_conformal/gpu_utilization_summary.json`
