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
