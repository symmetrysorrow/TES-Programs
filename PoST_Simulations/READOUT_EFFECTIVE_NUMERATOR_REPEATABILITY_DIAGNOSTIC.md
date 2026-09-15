# Effective numerator order-2 repeatability diagnostic

## Purpose

The effective-numerator complexity ladder showed that the 1--40 kHz repeat data can be represented well by a four-parameter order-2 model:

```text
P(x) = 1 + c2 x^2 + c4 x^4
x = f / 40 kHz
```

with a second-order pole denominator.

This diagnostic asks whether that reduced shape is stable within the 2024-12-05 run and whether the 2024-12-06 same-condition reference lies outside that within-run variation.

## Equal-record comparison

The 2024-12-06 reference accepted-record count is used as the block size.

With the current data:

```text
reference: 345 accepted records
repeat:    1731 accepted records
blocks:    5 x 345 records
remainder: 6 records
```

The final incomplete remainder is reported but excluded from block fits.

This makes the cross-day reference fit and every within-repeat block fit use the same number of records.

## Model

Every dataset is fit independently with order 2 only:

```text
transfer magnitude = sqrt(P(x) / D(f))

P(x) = (1 + u x^2)^2 + v^2 x^2
     = 1 + c2 x^2 + c4 x^4
```

The comparison parameters are:

- `pole_Hz`
- `pole_Q`
- `c2`
- `c4`

The latent optimization coordinates `u` and `v` are retained only for reproducibility. They are not used as physical or repeatability parameters.

## White-floor policy

The reference uses its detector-snapshot production post-filter white floor with scale 1.

The repeat full run and all repeat blocks use the same previously profiled repeat white-floor scale from:

```text
PoST_Simulations/cases/tagawa_20241205_r1ch12_215mK_1400uA_gain5_repeat/lowmid_free_transfer_snapshot.json
```

The white floor is not refit per block.

This keeps block-to-block changes from being hidden by a block-specific white nuisance amplitude.

## Fit and holdout

Every independent order-2 fit uses only:

```text
1 kHz <= f < 40 kHz
```

The strict holdout is:

```text
40 kHz <= f <= 200 kHz
```

No holdout sample is passed to the optimizer.

Each output fit contains fit-region, holdout, and full 1--200 kHz metrics.

## Parameter repeatability

For each canonical parameter the output reports:

- reference value;
- repeat-full value;
- all five repeat-block values;
- block median, mean, sample standard deviation, minimum, and maximum;
- whether the reference lies outside the repeat-block range;
- reference/repeat-full ratio.

For positive parameters it additionally reports log10-space block scatter and:

```text
abs(log10(reference) - median(log10(blocks)))
------------------------------------------------
sample std(log10(blocks))
```

This is a descriptive shift-to-scatter ratio, not a formal Gaussian significance.

The default screen marks a parameter shift when this ratio is at least 2.

## Transfer-shape repeatability

Parameter covariance can make individual coefficients move while the transfer shape remains nearly unchanged.

The diagnostic therefore also compares the normalized order-2 transfer shapes directly.

For each repeat block it computes the dB RMS difference from the repeat-full fitted transfer over:

- 1--40 kHz fit region;
- 40--200 kHz holdout region.

It then compares the reference-vs-repeat-full transfer-shape RMS against the largest repeat-block-vs-repeat-full RMS.

Important fields:

```text
transfer_shape_repeatability.reference_vs_repeat_full_fit_region
transfer_shape_repeatability.repeat_block_vs_repeat_full_fit_region
transfer_shape_repeatability.max_repeat_block_fit_region_rms_dB
transfer_shape_repeatability.day_fit_region_rms_over_block_max_ratio
```

The transfer-shape comparison is the preferred repeatability diagnostic because it is less sensitive than individual parameter values to correlated reparameterization.

## Interpretation flags

### `day_order2_transfer_shape_exceeds_within_repeat_block_variation`

True when the reference-vs-repeat-full 1--40 kHz transfer-shape RMS exceeds the largest same-run block-vs-full RMS.

### `at_least_two_canonical_parameters_shift_beyond_block_sigma_screen`

True when at least two of `pole_Hz`, `pole_Q`, `c2`, and `c4` exceed the configured 2-sigma descriptive log-scatter screen.

### `day_specific_order2_shape_shift_screen_passes`

True only when both the transfer-shape and multiple-parameter screens pass.

This is still not a physical electronics-drift identification.

The output always keeps:

```text
physical_electronics_drift_identified = false
```

## Git-tracked inputs

Configuration:

```text
PoST_Simulations/config/readout_effective_numerator_repeatability_config.json
```

Cross-dataset manifest:

```text
PoST_Simulations/config/shared_readout_cross_dataset_manifest.json
```

Reference transfer:

```text
PoST_Simulations/config/shared_readout_reference_transfer.json
```

Repeat white-floor/free-fit snapshot:

```text
PoST_Simulations/cases/tagawa_20241205_r1ch12_215mK_1400uA_gain5_repeat/lowmid_free_transfer_snapshot.json
```

Raw CH0 acquisition records remain external.

## Run

From the repository root:

```powershell
python PoST_Simulations/subScript/readout_effective_numerator_repeatability_diagnostic.py
```

Default output:

```text
PoST_Simulations/.noise_optimization_work_rsh_sweep/readout_effective_numerator_repeatability_diagnostic.json
```

## Next interpretation

If the reference transfer-shape difference is larger than every repeat block fluctuation and multiple canonical coefficients shift beyond the within-repeat scatter, the reduced order-2 shape shows a reproducible day-specific change.

If parameters move but the transfer shape remains inside the block envelope, the apparent coefficient drift is mainly parameter covariance.

If both parameters and shape vary strongly among the five repeat blocks, the dominant issue is within-run nonstationarity rather than a clean adjacent-day shift.
