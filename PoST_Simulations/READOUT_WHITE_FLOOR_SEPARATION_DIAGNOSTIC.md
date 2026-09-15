# Readout white-floor separation diagnostic

## Purpose

The adjacent-day repeatability work showed that the fixed reference transfer strongly improves the 2024-12-05 repeat, but a local transfer refit improves it further. The model-free day-to-day ASD comparison then showed that the fitted transfer-ratio shape does not reproduce the full measured day-to-day drift, especially above about 100 kHz.

This diagnostic separates high-frequency post-filter white-floor leverage from pre-ADC transfer leverage.

## Separation order

The diagnostic is deliberately ordered:

```text
Stage 0
shared reference transfer
+ inherited production post-filter white floor

Stage 1
shared reference transfer held completely fixed
+ profile only the post-filter white ASD amplitude

Stage 2
best Stage-1 white ASD held fixed
+ refit only the five infinity-pole hybrid transfer parameters
```

This order prevents the transfer fit from immediately using high-frequency transfer curvature to compensate for a different white floor.

## Git-tracked inputs

Standard execution requires no user-created JSON.

Configuration:

```text
PoST_Simulations/config/readout_white_floor_separation_config.json
```

Shared-transfer reference:

```text
PoST_Simulations/config/shared_readout_reference_transfer.json
```

Previous repeat-local transfer snapshot:

```text
PoST_Simulations/cases/tagawa_20241205_r1ch12_215mK_1400uA_gain5_repeat/local_readout_transfer_snapshot.json
```

The repeat detector snapshot and raw-record selection spec are obtained through the tracked shared-readout manifest.

Raw CH0 records remain on the acquisition drive.

## White-floor profile

The Stage-0 white ASD is the repeat detector snapshot's production post-filter white floor after the normal `post_filter_white_fraction` conversion.

Stage 1 multiplies that ASD by a non-negative scale while keeping detector parameters and the shared transfer fixed.

The tracked default profile covers:

```text
exact zero
0.01 x ... 30 x baseline on a logarithmic grid
baseline 1 x explicitly
+ bounded scalar optimization inside the same positive interval
```

The best candidate is selected by the same production shape score used by the surrounding diagnostics.

Important fields:

```text
white_floor_profile.baseline_white_asd_A_rtHz
white_floor_profile.best.white_scale
white_floor_profile.best.white_asd_A_rtHz
white_floor_profile.best.score_ratio_to_scale_one
stage1_shared_transfer_profiled_white.score_ratio_to_stage0
```

Absolute CH0 calibration remains unresolved, so the fitted white ASD is a diagnostic nuisance amplitude, not an identified physical source.

## Transfer refit after profiling white

Stage 2 freezes the best Stage-1 white ASD and refits only:

- second-order pole frequency;
- second-order pole Q;
- second-order zero frequency;
- second-order zero Q;
- lead zero frequency.

The lead pole remains at the exact infinity limit.

The exact shared transfer is retained as a warm-start candidate, so the local Stage-2 result cannot be worse than the Stage-1 shared-transfer result except for tiny numerical tolerance.

Important fields:

```text
stage2_profiled_white_transfer_refit.shape_score
stage2_profiled_white_transfer_refit.shared_over_local_score_ratio
stage2_profiled_white_transfer_refit.shared_within_local_score_tolerance
```

## Did white-floor freedom remove the previous transfer drift?

The diagnostic compares the Stage-2 transfer parameters with:

1. the reference shared transfer;
2. the previously tracked 2024-12-05 local transfer.

It reports RMS and maximum log10 parameter-ratio drift before and after the white-floor profile.

Important fields:

```text
transfer_drift_comparison.previous_local_snapshot.drift_from_shared
transfer_drift_comparison.after_white_profile.drift_from_shared
transfer_drift_comparison.rms_log10_drift_after_over_before
transfer_drift_comparison.drift_reduced_after_white_profile
```

## Interpretation flags

### `white_floor_materially_improves_shared_fit`

True when profiling only the white floor reduces the shared-transfer score by at least 20 percent.

### `shared_transfer_close_to_local_after_white_profile`

True when the Stage-1 shared-transfer score is within the configured factor 1.25 of the Stage-2 local-transfer score.

### `transfer_parameter_drift_reduced_after_white_profile`

True when the Stage-2 RMS log10 parameter drift from the shared transfer is smaller than the previously tracked repeat-local drift.

### `supports_white_floor_absorbing_previous_transfer_drift`

True only when all three conditions above are satisfied.

This is a diagnostic screen, not a physical source-identification test.

The output always keeps:

```text
physical_white_source_identified = false
electronics_transfer_drift_identified = false
```

## Run

From the repository root:

```powershell
python PoST_Simulations/subScript/readout_white_floor_separation_diagnostic.py
```

Default output:

```text
PoST_Simulations/.noise_optimization_work_rsh_sweep/readout_white_floor_separation_diagnostic.json
```

## Interpretation

If Stage 1 removes most of the shared/local gap and Stage-2 parameters move back toward the reference shared values, the earlier local transfer drift was substantially compensating for high-frequency white-floor differences.

If Stage 1 improves only the 100--200 kHz tail but Stage 2 still requires approximately the same 3.75 kHz / 12.8 kHz / 64.6 kHz local transfer structure, then a genuine run-dependent low/mid-frequency dynamic remains after white-floor separation.

If white-floor profiling does little at all, the previous local-transfer improvement cannot be attributed mainly to a simple post-filter white nuisance amplitude.
