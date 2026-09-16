# Phase24 preconditioner-lagging isolation

Generated: `2026-09-15T17:08:50.974150+00:00`

## Controlled change

- BDF order: `2` (unchanged)
- `Phase24 HYPRE Reuse`: `True` (unchanged/enabled)
- `Phase24 Preconditioner Lagging`: `disabled` (only intended change)
- Window: approximately `1 us` after the pulse

## Run result

- Process exit code: `0`
- `ALL DONE` marker: `True`
- Abort/allocation markers: `0`
- HYPRE solves: `45`; zero-iteration solves: `22`
- Relative-change samples: `45`; zero values: `16`

`ALL DONE` is a process-completion marker, not proof that every nonlinear solve converged. Review the solver log and the zero-iteration count together.

Comparison: `D:\Github\TES-Programs\Elmer-Projects\artifacts\comparison\phase24_bdf2_reuse_on_lagging_disabled_1us` (exit `0`)

Solver log: `D:\Github\TES-Programs\Elmer-Projects\results\case_phase24_short_bdf2_reuse_on_lagging_disabled_1us\solver.log`
Generated project: `D:\Github\TES-Programs\Elmer-Projects\artifacts\phase24_lagging_isolation_1us\phase24_lagging_isolation.json`
