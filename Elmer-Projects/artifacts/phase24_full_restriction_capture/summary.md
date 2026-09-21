# Phase24 full restriction capture

The diagnostic was enabled only for the pulse-off, BDF1, one-step mortar-on
case.  It captures `CollectionMatrix`, `CollectionVector`, and
`CollectionSolution` inside `SolveWithLinearRestriction`, immediately around
`SolveLinearSystem`.  Physics, mesh, material, timestep, and solver settings
were unchanged.

Runtime provenance for nl1–nl3 was identical:

| quantity | value |
| --- | ---: |
| physical time | 0.020018 s |
| dt / BDF order | 18 us / 1 |
| solver / direct method | direct / MUMPS |
| total / primal / constraint rows | 96955 / 96769 / 186 |
| CollectionMatrix nnz | 1329477 |
| constraint matrix rows | 186 |
| eliminate / penalty / restriction active | 0 / 0 / 1 |
| D block | exactly zero in saved CRS data |

Independent SuperLU results are in `full_system_analysis.json`.  Per
iteration, the direct full-system residual was `6.75–7.43e-11` relative and
the native full-system residual was `1.35–1.37e-10` relative.  Native versus
direct differences were at most `1.84e-7 K` in the primal vector and
`5.20e-7` in the full vector.  The block residuals were approximately
`7e-11` primal and `5–11e-11` constraint relative to the matrix/solution
scale; the constraint RHS is zero, so RHS-relative constraint residual is not
defined.

For nl1, the TES averages were:

| solution | TES average (K) |
| --- | ---: |
| x_before | 0.168569130247802 |
| primal-only K^-1 b | 0.166569454054102 |
| full direct primal | 0.167800479716859 |
| full native primal | 0.167800482506837 |
| HeatSolve outer_after | 0.167800482680851 |

The full native primal component agrees with `HeatSolve outer_after` within
`5.0e-8 K` maximum over primal DOFs.  Its TES-average difference is below
`6.3e-9 K` for nl1–nl3.

The mortar correction is the key result: for nl1,

`full_direct_primal - primal_only_direct` = `+1.231025663 mK` in the TES
weighted average, reproducing the previously unexplained `+1.231029 mK`
offset.  The direct `B^T lambda` norm was `1.2551e-11`, with maximum entry
`1.1851e-12`; its TES-weighted average contribution is `-2.58e-13` in the
assembled equation units.  The explicit D block has zero nnz and zero norm.

The before/after raw matrix and RHS hashes are different, but only at print
precision: maximum absolute changes were `1.11e-16` in A and `8.47e-22` in b.
The saved matrix record count is `1194063`; the runtime CRS storage size is
`1329477` including structural zero entries.

Conclusion: this is Case A.  The prior primal-only `K x = b` comparison was
not the system solved by MUMPS.  The `-0.768643 mK` nl1 jump is required by
the full mortar-constrained system, not evidence of a MUMPS result-vector or
post-solve transfer failure.  A small (~50 nK max) outer variable transfer
difference remains and should be checked next in the exact post-
`DefaultSolve` copy/reconstruction path, but it is not the source of the
millikelvin-scale discrepancy.

The raw `full_*.dat` files are retained locally under each `ts0001_nl000{1,2,3}`
directory; they remain ignored as large generated data.  `metadata.json`
contains the provenance and SHA-256 values for each capture.
