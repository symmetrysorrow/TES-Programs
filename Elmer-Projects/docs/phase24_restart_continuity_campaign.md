# Phase24 restart-continuity campaign

This diagnostic localizes the nonconforming 32-layer + mortar discontinuity
between the Gate3 steady state and the first pre-pulse transient state.

The target symptom is:

- Gate3 steady current: approximately 143.567589 uA
- MUMPS 40 us transient pre-pulse baseline: approximately 165.7 uA
- the first transient series point is already inconsistent before the pulse
  response can be evaluated.

The campaign therefore does not run another 40 us or 100 us production trace.
It audits restart provenance and runs short pulse-OFF controls.

## Run

From \`Elmer-Projects\`:

\`\`\`powershell
python scripts/support/run_phase24_restart_continuity_campaign.py
\`\`\`

If the current Gate3 / Gate4-5 case definitions are in a local generated
project that is not committed to main, pass it explicitly:

\`\`\`powershell
python scripts/support/run_phase24_restart_continuity_campaign.py \`
  --project path\\to\\project.json
\`\`\`

Useful isolation switches:

\`\`\`powershell
# Inspect existing artifacts only.
python scripts/support/run_phase24_restart_continuity_campaign.py --audit-only

# Generate the project and commands without starting Elmer.
python scripts/support/run_phase24_restart_continuity_campaign.py --dry-run

# Test the T0/fallback circuit initialization path too.
python scripts/support/run_phase24_restart_continuity_campaign.py --include-state-fallback

# Add no-mortar matrix-capture + observable hold controls. Diagnostic only on a nonconforming mesh.
python scripts/support/run_phase24_restart_continuity_campaign.py --include-no-mortar
\`\`\`

## Default short variants

The source transient case is copied and pulse energy is set to zero. Every
variant has a unique case name and does not overwrite production results.

| variant | backend | BDF | steps | mortar |
|---|---|---:|---:|---|
| current-state first step | MUMPS | 1 | 1 | on |
| Gate3-state first step | MUMPS | 1 | 1 | on |
| Gate3-state hold | MUMPS | 1 | 5 | on |
| Gate3-state production-integrator hold | MUMPS | production | 5 | on |
| Gate3-state backend comparison | HYPRE | 1 | 1 | on |
| Gate3-state no-mortar matrix capture (`--include-no-mortar`) | MUMPS | 1 | 1 | off |
| Gate3-state no-mortar observable hold (`--include-no-mortar`) | MUMPS | 1 | 5 | off |

## TES state-file safety

\`tes_transient_heat_source.f90\` can checkpoint \`TES State File\` after an
accepted transient timestep. A diagnostic must therefore never point multiple
variants at the original Gate3 state file, or even at one shared copied file.

The campaign creates a separate working state copy for every variant under:

\`\`\`text
artifacts/phase24_restart_continuity_campaign/state_snapshots/
\`\`\`

It records the five-value state before launching Elmer, because the working
copy can legitimately change after an accepted step.

A steady and transient case that use the same original \`TES State File\` path
are flagged as mutable restart provenance. A previous transient can otherwise
change the seed seen by a later run.

## Artifacts

\`\`\`text
artifacts/phase24_restart_continuity_campaign/
  summary.md
  summary.json
  provenance.json
  input_diff.json
  restart_audit.json
  state_file_audit.json
  case_matrix.csv
  first_step_metrics.csv
  field_continuity.json
  linear_system_comparison.json
  old_good_route_diff.json
  phase24_restart_continuity_campaign.json
  state_snapshots/
  dumps/
\`\`\`

TES series rows are committed when the following timestep begins. Therefore a one-step diagnostic can finish successfully with `row_count=0`; it is retained for matrix capture, while the hold5 cases are used for accepted-step BDF/mortar diagnosis.

The summary classifies the earliest captured divergence as:

1. TES state file before the solver,
2. first nonlinear / iteration record,
3. first accepted timestep,
4. not observed in the captured checkpoints.

The first question is continuity, not COMSOL parity:

\`\`\`text
Gate3 steady
  -> loaded thermal restart + circuit seed
  -> first nonlinear solve
  -> first accepted pulse-OFF timestep
\`\`\`

Only after this chain holds near the Gate3 operating point should the full
COMSOL transient waveform be used as the solver/physics parity gate.


## First-step restart residual decomposition

For matrix-capture variants the campaign now evaluates the saved first-solve
candidate against the assembled system:

```text
r = b - A x_saved
```

Rows with a nonzero diagonal are classified as primal thermal rows. Rows
without a nonzero diagonal are classified as mortar constraint rows. This
matches the explicit saddle structure currently observed on the nonconforming
32-layer path.

The Elmer SaveLinearSystem solution file is treated deliberately as an
x0/restart **candidate**, not as proof of the exact native pre-solve vector.
Existing Phase24 captures save only the primal temperature vector while the
mortar multiplier tail can be absent. Missing solution entries are therefore
extended with zero and this caveat is recorded in the artifact.

The output is:

```text
artifacts/phase24_restart_continuity_campaign/
  first_step_restart_residual.json
```

For each captured variant it reports full, primal, and constraint residual
L2/max norms plus a backward-error style normalization:

```text
||r|| / (||A||_F ||x|| + ||b||)
```

Interpretation:

- constraint block disproportionately large: inspect restart DOF
  reconstruction / mortar initialization;
- primal block large while constraint block is small: inspect steady-vs-
  transient heat assembly, source, BC, or material evaluation;
- both blocks small but the accepted TES state still jumps: inspect nonlinear
  electrothermal update and solution-transfer/bookkeeping paths.

The one-step mortar and no-mortar variants retain matrix/solution capture for
this diagnostic; the corresponding hold5 variants remain the source for
accepted-step physical continuity.


## Independent first-step direct solve and block comparison

The campaign now performs a best-effort independent sparse direct solve of the
saved first transient system when NumPy/SciPy are available:

```text
A x_direct = b
delta = x_direct - x_saved
```

This is intentionally separate from the Elmer/MUMPS solve. It answers whether
the *saved first linear system itself* requires a material correction from the
saved restart/x0 candidate.

Outputs:

```text
first_step_matrix_blocks.json
first_step_direct_solve.json
primal_system_comparison.json
```

The matrix is partitioned algebraically as:

```text
A = [ K   Bt ]
    [ B    D ]
```

using nonzero-diagonal rows as primal thermal rows and zero-diagonal rows as
constraint rows. The block artifact reports row counts, record counts,
Frobenius norms, maximum coefficients, B-vs-Bt-transpose mismatch, and primal
versus constraint RHS norms.

When `--include-no-mortar` is used, `primal_system_comparison.json`
compares only the primal K block and primal RHS of the mortar and no-mortar
one-step systems. If K and the primal RHS are unchanged, the first-step
difference is isolated to the interface constraint coupling rather than the
underlying thermal operator.

For the independent direct solution, the key quantity is:

```text
primal.delta_max_abs_mK_if_temperature
```

Interpretation:

- material direct-solve temperature correction: the assembled first transient
  linear system itself moves the saved restart/x0 state;
- tiny direct-solve correction but large accepted TES jump: focus on later
  nonlinear electrothermal coupling or solution-transfer/bookkeeping;
- identical primal K/RHS with mortar/no-mortar: focus on B/Bt/D coupling;
- changed primal K or RHS: Apply Mortar BCs is changing more than an appended
  constraint block on this path.

The direct solve is optional at runtime. If SciPy is unavailable or the sparse
factorization fails, the campaign records that status and still produces the
matrix-block and primal-system diagnostics.


## Nonlinear linear-system sequence

The restart campaign now inspects all continuously numbered SaveLinearSystem
files generated by the one-step Gate3-state MUMPS diagnostic instead of using
only the first dump. The first three saved systems are compared in dispatch
order:

```text
solve1 -> solve2
solve2 -> solve3
```

The output is:

```text
artifacts/phase24_restart_continuity_campaign/
  nonlinear_linear_system_sequence.json
```

For each adjacent pair it reports full matrix/RHS differences and, when the
primal/constraint partition is unchanged, block-resolved differences:

```text
A:
  K
  Bt
  B
  D

b:
  primal
  constraint
```

Each block includes L2 difference, maximum absolute difference, number of
changed records, and relative L2 difference with respect to the earlier
system.

This diagnostic is motivated by the observed first-timestep sequence:

```text
nonlinear iter 1 -> 2:
  TES temperature remains continuous
  relaxed Joule power changes by only a few ppm

nonlinear iter 2 -> 3:
  TES temperature drops by about 0.77 mK
```

Therefore the main question is whether the saved thermal system changes before
that temperature jump. Interpretation:

- RHS change dominates: inspect transient-history RHS, body-force insertion,
  and RHS-only/reuse assembly;
- operator change dominates: inspect temperature-dependent material assembly,
  mortar reconstruction, and matrix reuse;
- K is stable but B/Bt changes: interface constraint reconstruction is the
  leading target;
- matrix blocks are stable but primal RHS changes: transient mass/history or
  body-source RHS is the leading target.

The dump ordinal is intentionally described as SaveLinearSystem dispatch
order. On this single-step direct-MUMPS diagnostic it is expected to track the
nonlinear linear solves, but the artifact does not pretend the filename itself
contains an authoritative Elmer nonlinear-iteration ID.


## Direct sensitivity of the nonlinear solve transition

The campaign now solves four combinations for each adjacent saved-system pair:

```text
x11 = A1^-1 b1
x12 = A1^-1 b2
x21 = A2^-1 b1
x22 = A2^-1 b2
```

This separates the observed saved-system transition into:

```text
full change       = x22 - x11
RHS-only effect   = x12 - x11
operator-only     = x21 - x11
RHS effect on A2  = x22 - x21
A effect on b2    = x22 - x12
interaction       = x22 - x12 - x21 + x11
```

The artifact is:

```text
artifacts/phase24_restart_continuity_campaign/
  nonlinear_transition_direct_sensitivity.json
```

For scalar HeatSolve primal rows, all effects include both L2 and maximum
temperature-equivalent changes in mK. The headline reports the solve1->solve2
full, RHS-only, and operator-only maximum primal changes.

The decomposition uses two sparse LU factorizations, one for A1 and one for A2,
then reuses them for the crossed right-hand sides. It is therefore much more
diagnostic than comparing raw relative A/b norms alone: a numerically small
RHS perturbation can still produce a large temperature correction in an
ill-conditioned transient saddle system.

Interpretation:

- RHS-only temperature effect nearly equals the full transition while the
  operator-only effect is tiny: focus on transient-history/body-force/RHS
  assembly and reuse;
- operator-only effect nearly equals the full transition: focus on material,
  mortar, and matrix assembly/reuse;
- both are material: inspect the reported interaction and the block-resolved
  A/b differences before assigning one subsystem as the cause.

As with the earlier independent direct solve, this diagnostic is best-effort.
If SciPy is unavailable or sparse LU fails, the campaign records the reason
and continues with the non-direct diagnostics.


## Matrix dimension and zero-RHS mortar rows

SaveLinearSystem is configured with:

```text
Linear System Save Skip Zeros = True
```

For an explicit mortar saddle system,

```text
A = [ K   Bt ]
    [ B    0 ]

b = [ f ]
    [ 0 ]
```

the multiplier rows can therefore be absent from the saved RHS even though
they are present in the assembled matrix. Using `max(rhs_index)` as the
system dimension truncates those rows and can make two equal-size nonlinear
systems appear to have different dimensions.

The campaign now derives the dimension from the saved matrix A:

```text
n = max(max(A row index), max(A column index))
```

and zero-extends omitted RHS and solution entries to that dimension. This rule
is used by restart residual decomposition, K/B/Bt/D decomposition, independent
direct solve, primal-system comparison, nonlinear A/b comparison, and the
crossed direct sensitivity solve.

The campaign also discovers `*_sizes.dat` files and records their integer
contents as provenance metadata, but does not assume a build-specific semantic
interpretation. A-derived index extent is the authoritative diagnostic
dimension.

Relevant output fields now include:

```text
dimension_source
dimension.rows
rhs_saved_records
rhs_implicit_zero_entries
sizes_metadata
```

This also means earlier constraint-row counts derived from RHS extent should
be treated as provisional. Re-running the campaign after this change is
required before interpreting the multiplier-row count or a previous
"different dimensions" direct-sensitivity failure.
