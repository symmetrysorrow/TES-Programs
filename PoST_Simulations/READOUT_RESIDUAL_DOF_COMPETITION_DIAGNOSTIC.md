# Residual readout degree-of-freedom competition

## Purpose

The stability-aware R_TES profile showed that broad R_TES freedom plus the tested detector-state nuisance parameters still does not reach the 2024-12-05 repeat-local order-2 readout fit.

This diagnostic asks the next minimal question:

> After giving the detector/DC operating-point confounds their best tested freedom, which effective readout-shape degree of freedom is still required?

## Detector/DC confound baseline

The diagnostic starts from the Git-tracked best stability-aware profile point:

```text
R/R0 = 1.25
R_TES = 21.93979072046886 mOhm
```

with the tracked detector nuisance solution near:

```text
alpha  = 181.616
beta   = 0.00517
C_tes  = 3.4418e-13 J/K
L      = 12.29997 nH
T_bath = inherited 0.215701956 K
```

The fixed-reference-readout score at that baseline is:

```text
shape score = 2.1461327801482766e-4
RMS residual = 0.2544565 dB
```

The repeat-local order-2 comparator remains:

```text
shape score = 3.309316879681855e-5
RMS residual = 0.1104662 dB
```

The baseline is stored in:

```text
PoST_Simulations/cases/tagawa_20241205_r1ch12_215mK_1400uA_gain5_repeat/stability_aware_rtes_profile_snapshot.json
```

## Important fairness rule

`R_TES` is fixed to the profiled 1.25 R0 value for the entire diagnostic.

However, every residual-readout family **re-optimizes**:

```text
alpha
beta
C_tes
L
T_bath
```

This avoids overestimating readout freedom by freezing detector nuisance values that might redistribute once a readout coordinate is released.

The post-filter white ASD remains fixed to the tracked 2024-12-05 profiled value.

The detector nuisance bounds are also stored explicitly in the config, including the alpha fraction/max, beta range, T_bath half-width, and the use of production C_tes/L bounds.

## Readout parameterization

The readout is kept in the reduced order-2 effective form:

```text
P(x) = 1 + c2*x^2 + c4*x^4
x = f / 40 kHz
```

with a second-order pole denominator.

The canonical readout coordinates are:

```text
pole_Hz
pole_Q
c2
c4
```

`c2` and `c4` are effective magnitude-shape coefficients. They are not physical zero/component identifications.

## Main nested ladder

The primary comparison is strictly nested:

### 0. `fixed_reference_readout`

No readout parameter is free. Only the five detector nuisance parameters are refit.

### 1. `pole_frequency_only`

Release:

```text
pole_Hz
```

### 2. `pole_section`

Release:

```text
pole_Hz
pole_Q
```

### 3. `pole_section_plus_c2`

Release:

```text
pole_Hz
pole_Q
c2
```

### 4. `full_order2`

Release:

```text
pole_Hz
pole_Q
c2
c4
```

Thus the main ladder asks whether pole frequency alone is sufficient, whether Q adds necessary leverage, whether the primary numerator curvature `c2` is still required, and whether any residual need for `c4` remains.

## Branch diagnostics

Two additional non-ladder branches isolate numerator leverage:

### `c2_only`

Release only `c2` while keeping the entire reference pole section and `c4` fixed.

### `pole_frequency_plus_c2`

Release:

```text
pole_Hz
c2
```

while keeping `pole_Q` and `c4` at the reference values.

This branch tests whether the main residual information can be represented by the robust pole-frequency scale plus the primary numerator curvature without requiring Q freedom.

## Optimizer

Every family uses differential evolution followed by least-squares refinement when stable.

Warm-start candidates include:

- the tracked detector/DC baseline with reference readout;
- the tracked detector/DC baseline with repeat-local readout values projected onto the free coordinates;
- previously completed residual-readout family solutions projected onto the current free coordinates.

Each detector trial is rejected unless the TES operating point is valid and stable.

## Fit and holdout

Optimization uses only:

```text
1 kHz <= f < 40 kHz
```

The strict holdout is:

```text
40 kHz <= f <= 200 kHz
```

No holdout point enters the optimizer.

Every family reports fit-region, holdout, and full 1--200 kHz metrics.

## Regression anchor

`fixed_reference_readout` must reproduce the tracked R=1.25 detector/DC baseline.

The diagnostic aborts if its score is more than 0.05 percent worse than the tracked baseline.

## Near-local screen

A residual-readout family reaches the repeat-local benchmark only if both conditions hold:

```text
family score / repeat-local score <= 5
family RMS - repeat-local RMS <= 0.05 dB
```

This is the same practical screen used by the preceding confound diagnostics.

## Incremental output

For the main nested ladder, each step reports:

```text
shape_score_ratio_new_over_old
rms_improvement_dB
holdout_score_ratio_new_over_old
material_score_improvement
```

The default material-improvement screen is a score ratio <= 0.8.

## Key fields

```text
residual_readout_fits

main_nested_ladder.incremental_gain
main_nested_ladder.minimal_family_near_repeat_local

branch_diagnostics.c2_only
branch_diagnostics.pole_frequency_plus_c2
branch_diagnostics.minimal_any_family_near_repeat_local
```

Important flags:

```text
interpretation_flags.pole_frequency_only_sufficient
interpretation_flags.pole_section_sufficient
interpretation_flags.c2_only_sufficient
interpretation_flags.pole_frequency_plus_c2_sufficient
interpretation_flags.pole_section_plus_c2_sufficient
interpretation_flags.full_order2_reaches_repeat_local
interpretation_flags.c4_required_after_pole_section_plus_c2
```

## Interpretation

If `pole_frequency_only` reaches the comparator, the residual day-specific information after detector/DC freedom is concentrated mainly in the pole-frequency scale.

If `c2_only` succeeds while pole-only does not, the primary numerator curvature is the more independent residual shape lever.

If `pole_frequency_plus_c2` succeeds but `pole_section` does not, the data favor pole scale plus numerator curvature over an independently shifted pole Q.

If `pole_section_plus_c2` succeeds and `full_order2` adds little, `c4` is not required after detector/DC confounds are allowed.

If only `full_order2` succeeds, the residual low/mid shape still requires the complete reduced order-2 effective readout family.

None of these outcomes identifies a physical electronics component. They identify only the minimum effective magnitude-shape freedom that remains after the tested detector/DC nuisance set.

## Git-tracked inputs

Configuration:

```text
PoST_Simulations/config/readout_residual_dof_competition_config.json
```

Baseline snapshot:

```text
PoST_Simulations/cases/tagawa_20241205_r1ch12_215mK_1400uA_gain5_repeat/stability_aware_rtes_profile_snapshot.json
```

Cross-dataset manifest:

```text
PoST_Simulations/config/shared_readout_cross_dataset_manifest.json
```

Raw CH0 records remain external.

## Run

From the repository root:

```powershell
python PoST_Simulations/subScript/readout_residual_dof_competition_diagnostic.py
```

Default output:

```text
PoST_Simulations/.noise_optimization_work_rsh_sweep/readout_residual_dof_competition_diagnostic.json
```
