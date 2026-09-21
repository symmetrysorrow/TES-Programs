# Phase24 restart-continuity campaign

Generated: `2026-09-20T15:57:21.844932+00:00`

## Headline

- First observed 143.57 -> abnormal-current divergence: **first accepted timestep**
- Restart File audit: **present**
- TES State File audit: **matches current gate at file level; steady/transient state hashes are identical**
- Mortar initialization implicated: **mortar-dependent first-step divergence observed**
- Time-integration initialization implicated: **no material first-step BDF effect observed**
- MUMPS/HYPRE backend difference implicated: **not isolated**
- Strongest current explanation: **only one SaveLinearSystem system was captured for the mortar case; textual sibling variants were previously misclassified as numbered dumps. The first linear system remains near the restart state, while the later thermal jump still requires a true per-nonlinear-iteration outer-system capture**
- Series observability: TES accepted-step series is written when the next timestep begins; one-step variants may legitimately have row_count=0. Hold5 variants are used for BDF and mortar first-step diagnosis.
- Restart/x0 residual: **captured: constraint rows=186, constraint backward error=2.797506967746546e-10, primal backward error=2.554113637668745e-10**
- Matrix dimension audit: **A-derived rows=96955, saved RHS records=96955, implicit-zero RHS entries=0, constraint rows=186**
- Independent direct first solve: **captured: primal max |x_direct-x_saved|=0.0001676485824053664 mK, direct relative residual=7.416481296013727e-11**
- Mortar/no-mortar primal system: **same primal row set=True; K max diff=1.0946519945635895e-09; primal RHS max diff=1.8014562857879903e-13**
- Linear dump provenance: **dumps=1, same-shape outer candidates=1, heterogeneous=0, ignored sibling variants=1, log source=solver.log, log save events=1**
- Nonlinear A/b sequence: **not captured**
- Direct solve1->2 sensitivity: **not compared: fewer than two structurally homogeneous outer-candidate dumps**

## Evidence

- Reference Gate3 current: 143.567589 uA
- Source project: `artifacts\phase24_restart_continuity_campaign\phase24_restart_continuity_campaign.json`
- Existing transient case: `case_phase24_g45_s32m_40us_mumps_mortar_mumps`
- Existing steady case: `case_phase24_g3_s32m_s32m2_10f17b`

Detailed machine-readable outputs are beside this file: provenance.json, input_diff.json, restart_audit.json, state_file_audit.json, case_matrix.csv, first_step_metrics.csv, field_continuity.json, linear_system_comparison.json, first_step_restart_residual.json, first_step_matrix_blocks.json, first_step_direct_solve.json, primal_system_comparison.json, linear_dump_provenance.json, nonlinear_linear_system_sequence.json, nonlinear_transition_direct_sensitivity.json, and old_good_route_diff.json.

The campaign deliberately does not run a 40 us or 100 us production trace.
