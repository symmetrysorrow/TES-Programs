# Low/mid transfer fit with high-frequency holdout

## Purpose

The white-floor separation diagnostic showed that a higher post-filter white floor explains a substantial part of the repeat-case high-frequency mismatch, but a full-band local transfer fit still moved away from the shared reference transfer.

This diagnostic asks whether those remaining transfer shifts are actually demanded by the low/mid-frequency data.

The key rule is strict:

```text
fit:      1 kHz <= f < 40 kHz
holdout: 40 kHz <= f <= 200 kHz
```

No holdout point is passed to the optimizer.

## Git-tracked inputs

Configuration:

```text
PoST_Simulations/config/readout_lowmid_holdout_config.json
```

Reference shared transfer:

```text
PoST_Simulations/config/shared_readout_reference_transfer.json
```

White-profiled full-band repeat snapshot:

```text
PoST_Simulations/cases/tagawa_20241205_r1ch12_215mK_1400uA_gain5_repeat/white_profiled_full_band_transfer_snapshot.json
```

The snapshot freezes the preceding result:

- best white scale = 1.3635953106605805;
- best white ASD = 5.994841960426966e-11 A/rtHz;
- full-band local pole = 12.373 kHz;
- full-band local zero = 75.531 kHz;
- full-band local lead zero = 4.069 kHz.

These are diagnostic derived values, not hardware calibration constants.

## Fit procedure

1. Reconstruct the repeat pre-analysis target from the tracked accepted-mask spec.
2. Freeze detector parameters.
3. Freeze the post-filter white floor to the tracked profiled value.
4. Restrict model, target, and context arrays to 1--40 kHz.
5. Fit only the five infinity-pole hybrid transfer parameters.
6. Reconstruct the fitted model over the full 1--200 kHz band.
7. Evaluate 40--200 kHz only after the optimizer has finished.

The shared transfer is included as an exact warm-start candidate.

## Comparators

The output compares three transfers at the same fixed profiled white floor:

- reference shared transfer;
- preceding full-band local transfer;
- new 1--40 kHz-only transfer.

Metrics are reported separately for:

```text
low_mid_region_metrics
high_frequency_holdout_metrics
full_1_200k_metrics_for_low_mid_fit
```

## Parameter-return test

The main parameter-identifiability question is whether excluding 40--200 kHz causes the fitted parameters to move back toward the shared reference values.

For each transfer parameter the output reports absolute log10 distance from the shared value for both the full-band and low/mid-only fits.

Important fields:

```text
parameter_return_toward_shared.parameters
parameter_return_toward_shared.n_parameters_closer_to_shared
parameter_return_toward_shared.fraction_parameters_closer_to_shared
parameter_return_toward_shared.rms_log10_drift_low_mid_over_full_band
```

The zero-frequency check is highlighted separately:

```text
interpretation_flags.zero_frequency_returns_toward_shared_when_high_is_held_out
```

If the full-band zero near 75.5 kHz moves back toward the shared ~68.5 kHz value when high frequencies are excluded from fitting, that is evidence that the previous zero shift was driven by high-frequency leverage rather than by the low/mid data.

## Holdout generalization

After the 1--40 kHz fit is complete, its frozen parameters are evaluated on 40--200 kHz.

Important fields:

```text
high_frequency_holdout_metrics.shared_transfer
high_frequency_holdout_metrics.full_band_local_transfer
high_frequency_holdout_metrics.low_mid_fit_transfer
high_frequency_holdout_metrics.low_mid_over_shared_score_ratio
high_frequency_holdout_metrics.low_mid_over_full_band_score_ratio
```

The holdout comparison is descriptive. The full-band local fit has an expected advantage because it was allowed to use the holdout frequencies during fitting.

## Interpretation flags

### `low_mid_data_require_case_specific_transfer`

True when the low/mid-only local fit reduces the 1--40 kHz score by at least 20 percent relative to the shared transfer.

### `low_mid_fit_holdout_not_materially_worse_than_shared`

True when the low/mid-only transfer's 40--200 kHz holdout score is within the configured 1.25 factor of the shared transfer.

### `zero_frequency_returns_toward_shared_when_high_is_held_out`

True when the low/mid-only zero frequency is closer in log space to the shared reference zero than the preceding full-band local zero.

### `overall_transfer_drift_reduced_when_high_is_held_out`

True when the RMS log10 parameter drift from the shared transfer is smaller for the low/mid-only fit than for the full-band fit.

### `supports_full_band_zero_drift_being_high_frequency_nuisance`

True when both the zero frequency and the overall parameter vector return toward the shared transfer after high frequencies are excluded from fitting.

This remains diagnostic evidence only.

## Run

From the repository root:

```powershell
python PoST_Simulations/subScript/readout_lowmid_holdout_diagnostic.py
```

Default output:

```text
PoST_Simulations/.noise_optimization_work_rsh_sweep/readout_lowmid_holdout_diagnostic.json
```

## How to interpret the next result

If the low/mid-only fit still prefers roughly 4 kHz lead-zero and 12--13 kHz pole structure but the ~75 kHz zero returns toward ~68.5 kHz, the low/mid run-dependent dynamics are likely real while the higher-frequency zero shift was partly a nuisance response.

If all five parameters return close to the shared values, most of the previous local-transfer drift was likely induced by fitting the high-frequency region.

If the low/mid-only fit remains strongly displaced from the shared transfer even without any 40--200 kHz information, the run-dependent low/mid dynamics survive the white-floor and holdout tests.

If the low/mid-only fit performs badly on the high-frequency holdout, that is not by itself a failure: it means the high-frequency region contains additional structure not identified by the low/mid data.
