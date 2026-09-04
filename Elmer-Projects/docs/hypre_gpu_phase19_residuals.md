# Phase19 residual and acceptance definitions

All solver comparisons use the same assembled system and report the metrics
below separately.  A solver-reported relative residual is not interchangeable
with an absolute residual, a physical error, or a residual normalized by a
different right-hand-side norm.

| Metric | Definition | Availability rule |
|---|---|---|
| `solver_reported_relative_residual` | Value emitted by the iterative solver, using its own norm/scaling | Report exactly as emitted; do not relabel it |
| `original_system_absolute_residual` | `||A x - b||_2` for the original assembled system | Required for a full saddle vector |
| `original_system_relative_residual` | `||A x - b||_2 / ||b||_2` | Required for a full saddle vector |
| `constraint_absolute_residual` | `||B u + D lambda - g||_2` | Required for an explicit saddle vector |
| `constraint_relative_or_scaled_residual` | Constraint residual divided by `max(||Bu||, ||Dlambda||, ||g||)` | `null` when the scale is numerically zero |
| `relative_primal_error_vs_MUMPS` | `||u-u_MUMPS||_2 / ||u_MUMPS||_2` | Compare only primal entries |
| `TES temperature error` | Explicit physical temperature difference against the reference case | Report with its physical unit and norm/point definition |
| `TES current error` | Explicit current difference against the reference case | Report with its physical unit and norm/point definition |

For this case `g=0`, so `constraint_absolute_residual` is the primary
constraint gate.  A relative value against `||g||` is undefined and is not
reported.  The scaled value above is only a diagnostic.

The phrase “residual gate `1e-11`” refers to the original-system relative
residual unless a table explicitly says otherwise.  For example, an absolute
residual near `1e-16` can correspond to a relative residual near `3e-11` when
`||b||` is small; those are not the same metric.

The saved MUMPS reference contains `84,636` primal entries, not the full
`87,534`-entry saddle vector.  Its full saddle residual and multiplier
residual are therefore `not available`, not zero and not inferred.

## Current one-rank algebra table

These values are copied from the committed diagnostic artifacts in
`artifacts/hypre_phase19_schur/`.

| Solver | Full abs residual | Full rel residual | Constraint abs | Primal error vs MUMPS | TES temperature error | TES current error | Pass |
|---|---:|---:|---:|---:|---:|---:|---|
| Exact Schur oracle | `1.14e-16` | `3.18e-11` | `1.16e-24` | `2.21e-6` | not measured | not measured | algebraic reference |
| Diagonal Schur | `7.13e-6` | `1.99` | `4.98e-25` | `9.53e-1` | not measured | not measured | no |
| ILU proxy | `1.10e-5` | `3.06` | `6.58e-23` | `7.54e-1` | not measured | not measured | no |
| MGR best | not available in this artifact set | `2.12e-10` at iteration limit | not available | not accepted | not measured | not measured | no |

The exact oracle is a correctness reference, not a production strategy.  It
does not pass the `1e-11` relative gate in the saved double-precision report
because the absolute residual is close to machine precision while `||b||` is
small; it does satisfy the constraint algebra to the reported precision.

No production candidate is accepted until a full candidate reports all
available algebraic metrics and then passes the physical 7-step gate.
