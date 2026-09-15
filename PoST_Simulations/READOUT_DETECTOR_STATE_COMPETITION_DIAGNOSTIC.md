# Detector-state vs fixed-readout competition diagnostic

## Purpose

The order-2 repeatability diagnostic showed a reproducible 1--40 kHz effective-shape difference between the 2024-12-06 reference and the 2024-12-05 repeat, larger than the five equal-record repeat-block fluctuations.

That result still has an important confound: the 2024-12-05 detector snapshot is inherited from the 2024-12-06 operating-point fit.

This diagnostic therefore asks:

> Can a limited set of detector-state nuisance parameters absorb the 2024-12-05 low/mid shape while the 2024-12-06 order-2 readout shape is held fixed?

## Fixed quantities

The effective readout transfer is frozen to the Git-tracked 2024-12-06 order-2 solution:

```text
pole_Hz = 12300.466...
pole_Q  = 0.790743...
P(x) = 1 + c2*x^2 + c4*x^4
c2 = 42.0214
c4 = 34.2299
x = f / 40 kHz
```

The 2024-12-05 post-filter white ASD is frozen absolutely to:

```text
5.994841960426966e-11 A/rtHz
```

It is not recomputed as detector nuisance parameters move.

This prevents the detector-state optimizer from silently using white-floor freedom.

## Competing comparator

The diagnostic also reconstructs the tracked 2024-12-05 repeat-local order-2 solution with the inherited detector snapshot:

```text
pole_Hz = 11823.837...
pole_Q  = 0.657256...
c2 = 84.3812
c4 = 24.1775
```

This is the day-specific readout-shape comparator.

## Nested detector-state nuisance families

Three nested families are fit independently.

### 1. `transition_sensitivity`

Varies only:

```text
alpha
beta
```

### 2. `transition_plus_local_dynamics`

Varies:

```text
alpha
beta
C_tes
L
```

### 3. `transition_local_plus_bath`

Varies:

```text
alpha
beta
C_tes
L
T_bath
```

Everything else in the detector snapshot remains fixed, including `R_TES` and the thermal-link parameters.

## Bounds

The nuisance bounds reuse production conventions where available:

```text
alpha: 0.05 * inherited alpha ... 200
beta:  0 ... 12
C_tes: production C_TES_FIT_MIN ... C_TES_FIT_MAX
L:     production L_FIT_MIN ... L_FIT_MAX
T_bath: inherited value +/- 2 mK
```

`C_tes` and `L` are optimized logarithmically.

Every trial is screened with the TES operating-point validity/stability check before scoring.

## Fit and holdout

Detector nuisance optimization uses only:

```text
1 kHz <= f < 40 kHz
```

The strict holdout is:

```text
40 kHz <= f <= 200 kHz
```

No holdout point enters the optimizer.

Each nuisance family reports its fit-region, holdout, and full-band metrics.

## Main comparison

For every detector-state family the diagnostic reports:

```text
score_ratio_to_fixed_reference_inherited_detector
score_ratio_to_repeat_local_readout_comparator
rms_delta_to_repeat_local_readout_comparator_dB
```

The default `near_repeat_local_readout_comparator` screen requires both:

```text
detector_nuisance_score / repeat_local_readout_score <= 5
detector_nuisance_RMS - repeat_local_readout_RMS <= 0.05 dB
```

A detector nuisance family is considered materially useful relative to the inherited-detector/fixed-reference-readout baseline when:

```text
score ratio <= 0.5
```

## Boundary and operating-point diagnostics

Each nuisance fit reports:

- fitted values;
- ratio to inherited detector values;
- lower/upper boundary hits;
- TES operating-point validity/stability;
- current and Joule-power ratios relative to the inherited detector snapshot;
- optimizer instability rejects.

These fields are important because a mathematical nuisance fit that runs to boundaries or large operating-point changes is weaker evidence for a realistic day-state explanation.

## Interpretation flags

### `limited_detector_state_nuisance_can_compete_with_day_specific_readout_shape`

True only when the best tested detector-state family both materially improves the fixed-reference-readout baseline and reaches the near-repeat-local screen.

If true, the existing evidence does not distinguish a day-specific readout shape from detector-state drift within this limited nuisance family.

### `day_specific_readout_shape_still_required_within_tested_nuisance_set`

True when the tested detector-state nuisance families cannot reach the repeat-local readout comparator.

This means only that the tested nuisance set is insufficient.

It does **not** prove electronics drift.

## Important limitation

`R_TES` remains inherited/fixed because there is no independently linked 2024-12-05 IV operating-point fit in this diagnostic.

Therefore a negative detector-state result does not exclude an untested DC operating-point difference.

The output always keeps:

```text
physical_electronics_drift_identified = false
physical_detector_state_shift_identified = false
```

## Git-tracked inputs

Configuration:

```text
PoST_Simulations/config/readout_detector_state_competition_config.json
```

Order-2 day-comparison snapshot:

```text
PoST_Simulations/cases/tagawa_20241205_r1ch12_215mK_1400uA_gain5_repeat/order2_day_comparison_snapshot.json
```

Cross-dataset manifest:

```text
PoST_Simulations/config/shared_readout_cross_dataset_manifest.json
```

Raw CH0 records remain external.

## Run

From the repository root:

```powershell
python PoST_Simulations/subScript/readout_detector_state_competition_diagnostic.py
```

Default output:

```text
PoST_Simulations/.noise_optimization_work_rsh_sweep/readout_detector_state_competition_diagnostic.json
```

## Reading the result

If `alpha+beta` alone approaches the repeat-local comparator, the day-specific low/mid shape can plausibly be absorbed by transition-sensitivity state.

If only `+C_tes+L` works, the confound is specifically dynamical rather than a simple alpha/beta change.

If only `+T_bath` works, the result points toward operating-point/thermal-state sensitivity.

If no tested family approaches the repeat-local order-2 fit, the day-specific effective readout freedom remains necessary within this restricted detector-state model, and the next high-value test is an independently constrained 2024-12-05 operating point or IV-linked `R_TES` profile.
