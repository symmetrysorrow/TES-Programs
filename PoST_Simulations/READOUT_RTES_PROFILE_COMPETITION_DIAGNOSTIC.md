# R_TES profile vs detector-state/readout competition

## Purpose

The detector-state competition diagnostic showed that limited detector-state freedom can remove most of the 2024-12-05 versus 2024-12-06 low/mid mismatch, but it did not reach the repeat-local order-2 readout fit. The largest remaining untested operating-point confound is `R_TES`: the 2024-12-05 repeat currently inherits the 2024-12-06 IV-linked resistance.

This diagnostic profiles `R_TES` externally and refits the same nested detector-state nuisance families at every fixed resistance.

`R_TES` is never fitted simultaneously as a continuous optimizer coordinate.

## Fixed quantities

The following remain fixed throughout the profile:

- 2024-12-06 reference order-2 effective readout transfer;
- 2024-12-05 profiled absolute post-filter white ASD;
- all detector parameters not named by the selected nuisance family;
- exact accepted-record mask for the 2024-12-05 repeat.

Fit region:

```text
1 kHz <= f < 40 kHz
```

Strict holdout:

```text
40 kHz <= f <= 200 kHz
```

No holdout point is passed to the optimizer.

## R_TES profile

The inherited resistance is:

```text
R0 = 0.017551832576375086 ohm
```

The default fixed profile is:

```text
R/R0 = 0.75, 0.85, 0.90, 0.95, 1.00, 1.05, 1.10, 1.15, 1.25, 1.35, 1.50
```

The denser +/-15 percent region is intended to answer whether a modest DC operating-point shift is sufficient before interpreting more extreme profile points. The upper extension to 1.35 and 1.50 is included because the first profile had its best evaluated point at the former 1.25 upper edge.

## Stability-aware fixed-R semantics

A fixed-R point is **not** rejected merely because the inherited nuisance state at that resistance is unstable.

For every R point, the diagnostic records the inherited-state operating-point stability for provenance, then runs the nuisance optimizer. Every optimizer trial is screened by the TES validity/stability test.

Therefore:

```text
inherited nuisance state unstable
    !=
fixed R ruled out
```

A profile point is reported as:

```text
no_stable_solution_at_fixed_R
```

only when none of the configured nuisance-family searches produces a stable fitted candidate.

The output summarizes this explicitly under:

```text
R_TES_profile.stability_summary
interpretation_flags.low_R_stable_nuisance_solution_found
interpretation_flags.all_tested_R_points_have_stable_nuisance_solution
```

This change completes the low-R half of the profile that the original pre-screen could not evaluate.

## Nested nuisance refits

At each fixed R_TES point, the diagnostic independently refits:

### `transition_sensitivity`

```text
alpha, beta
```

### `transition_plus_local_dynamics`

```text
alpha, beta, C_tes, L
```

### `transition_local_plus_bath`

```text
alpha, beta, C_tes, L, T_bath
```

The nuisance bounds reuse the preceding detector-state competition conventions.

## R=1 regression anchor

The preceding detector-state result is stored in:

```text
PoST_Simulations/cases/tagawa_20241205_r1ch12_215mK_1400uA_gain5_repeat/detector_state_competition_snapshot.json
```

The tracked best R/R0=1 detector-state solution is supplied as an exact/warm-start candidate.

The diagnostic aborts if the R=1 profile best score is more than 0.05 percent worse than the tracked result.

This prevents the R-profile conclusion from being driven by a poorer optimizer basin than the already established R=1 fit.

## Comparators

The primary comparator is still the 2024-12-05 repeat-local order-2 readout fit with the inherited detector snapshot:

```text
shape score = 3.309316879681855e-05
RMS residual = 0.11046624584893718 dB
```

For every R/family row the output reports:

```text
score_ratio_to_repeat_local_readout_comparator
rms_delta_to_repeat_local_readout_comparator_dB
near_repeat_local_readout_comparator
```

The near-local screen is unchanged:

```text
score ratio <= 5
RMS penalty <= 0.05 dB
```

## Operating-point diagnostics

For each R point the output reports the fixed R value and the baseline TES operating point.

Every fitted nuisance family also reports:

- fitted nuisance values;
- boundary hits;
- current ratio to the profile-point baseline;
- Joule-power ratio to the profile-point baseline;
- optimizer instability rejects;
- 40--200 kHz holdout metrics;
- full 1--200 kHz metrics.

These fields are required for judging whether a mathematically improved R profile point is a credible operating-point alternative or another reduced-model boundary stress.

## Main output

Important fields are:

```text
R_TES_profile.rows
R_TES_profile.R1_regression_anchor
best_overall_profile_point
best_modest_R_shift_profile_point
```

`best_modest_R_shift_profile_point` is restricted to:

```text
0.90 <= R/R0 <= 1.10
```

## Interpretation flags

### `some_R_profile_point_reaches_near_repeat_local_readout_fit`

True when any tested fixed-R point plus detector nuisance reaches the established repeat-local order-2 readout screen.

### `modest_R_shift_reaches_near_repeat_local_readout_fit`

Same criterion, but restricted to +/-10 percent around the inherited R.

### `R_profile_materially_improves_over_R1_detector_nuisance`

True when the best R profile point reduces score to <=0.8 times the best R=1 detector-state nuisance score.

### `best_R_profile_point_is_at_tested_edge`

True when the best point is at the minimum or maximum tested R ratio. With the current profile that means 0.75 or 1.50 times R0. This is a warning that the selected profile range still did not bracket the optimum.

### `profiled_R_TES_can_compete_with_day_specific_readout_shape`

True when at least one fixed-R detector-state solution reaches the repeat-local readout comparator.

This would mean the current data cannot distinguish the day-specific effective readout shape from an unmeasured R_TES/DC operating-point shift plus the tested nuisance freedom.

### `day_specific_readout_shape_still_required_within_profiled_R_and_tested_nuisance_set`

True when no tested R point reaches the repeat-local comparator.

This strengthens the case that the effective day-specific shape is not merely the tested detector-state/DC operating-point confound, but it still does not identify physical electronics drift.

## Git-tracked inputs

Configuration:

```text
PoST_Simulations/config/readout_rtes_profile_competition_config.json
```

Order-2 day comparison:

```text
PoST_Simulations/cases/tagawa_20241205_r1ch12_215mK_1400uA_gain5_repeat/order2_day_comparison_snapshot.json
```

R=1 detector-state regression anchor:

```text
PoST_Simulations/cases/tagawa_20241205_r1ch12_215mK_1400uA_gain5_repeat/detector_state_competition_snapshot.json
```

Raw CH0 records remain external.

## Run

From the repository root:

```powershell
python PoST_Simulations/subScript/readout_rtes_profile_competition_diagnostic.py
```

Default output:

```text
PoST_Simulations/.noise_optimization_work_rsh_sweep/readout_rtes_profile_competition_diagnostic.json
```

## Reading the next result

If a modest R shift reaches the repeat-local comparator without severe nuisance-boundary stress, an independently linked 2024-12-05 IV/R_TES constraint becomes the highest-priority missing measurement.

If only an extreme edge R point works, the result is an existence proof rather than a plausible operating-point explanation.

If low-R inherited states are unstable but stable nuisance solutions are found after optimization, those fixed-R points must be judged by their fitted score, boundary stress, operating point, and holdout rather than by the inherited-state pre-screen.

If no R point works and the best nuisance solutions remain boundary-stressed, the day-specific effective readout freedom survives the major detector-state confounds tested so far.
