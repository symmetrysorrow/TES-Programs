# Model-free day-to-day readout drift diagnostic

## Purpose

The shared-readout cross-dataset diagnostic showed two things at once:

1. the fixed reference transfer strongly improves the adjacent-day 2024-12-05 repeat;
2. a local fit with the same infinity-pole hybrid topology improves the repeat further and shifts the fitted transfer parameters modestly.

The next question is whether that apparent transfer drift is visible directly in the measured data, without using the TES noise model.

This diagnostic therefore compares the **raw-record pre-analysis ASD shape** between the 2024-12-05 repeat and the 2024-12-06 reference.

The empirical day-to-day ratio is model-free.

## Git-tracked inputs

Standard execution needs no user-created JSON.

Configuration:

```text
PoST_Simulations/config/readout_day_to_day_drift_config.json
```

Cross-dataset manifest:

```text
PoST_Simulations/config/shared_readout_cross_dataset_manifest.json
```

Reference shared transfer:

```text
PoST_Simulations/config/shared_readout_reference_transfer.json
```

Tracked 2024-12-05 local-transfer snapshot:

```text
PoST_Simulations/cases/tagawa_20241205_r1ch12_215mK_1400uA_gain5_repeat/local_readout_transfer_snapshot.json
```

The local-transfer snapshot records the best local infinity-pole hybrid result from the preceding repeatability diagnostic. It is a derived diagnostic input, not a physical electronics identification.

Raw CH0 records remain on the acquisition drive.

## Empirical ASD construction

For each run:

```text
exact accepted CH0 records
-> per-record mean removal
-> Hann
-> power average
-> one-sided ASD
-> independent normalization at 1 kHz
```

No 10 kHz digital analysis filter is applied to the ASD used in the day-to-day comparison.

The empirical ratio is:

```text
normalized pre-analysis ASD(2024-12-05 repeat)
------------------------------------------------
normalized pre-analysis ASD(2024-12-06 reference)
```

No detector parameters, TES noise sources, thermal model, electrical model, or readout-transfer model are used to construct this ratio.

## Transfer-ratio overlay

Only after the empirical ratio has been formed, the diagnostic computes:

```text
normalized H_local_2024-12-05(f)
--------------------------------
normalized H_shared_2024-12-06(f)
```

where both transfer magnitudes use the same infinity-pole hybrid topology and are independently normalized at 1 kHz.

The empirical ratio and transfer ratio are compared through:

- 1--200 kHz dB-shape correlation;
- RMS of empirical minus transfer ratio in dB;
- broad-band mean/RMS residuals;
- anchor frequencies at 1, 5, 10, 20, 40, 70, 100, 150, and 200 kHz.

A good match means the fitted transfer drift follows an independently observed day-to-day spectral-shape change. It still does **not** prove that the electronics physically changed, because the detector operating point was not independently re-fit for the repeat.

## Block repeatability

The 2024-12-05 repeat contains substantially more accepted records than the reference run.

The diagnostic therefore uses the reference accepted-record count as the default block size. With the current data this is expected to be 345 records.

The 2024-12-05 accepted-record sequence is divided into deterministic, non-overlapping full blocks. Any final incomplete remainder is reported and excluded from the block envelope.

For every block, the pre-analysis ASD is independently normalized at 1 kHz and divided by the full 2024-12-05 normalized ASD.

The output reports:

- min / 16th / median / 84th / max block-ratio dB at anchor frequencies;
- maximum absolute block variation at each anchor;
- each block's 1--200 kHz RMS shape deviation from the full repeat ASD;
- maximum and median block RMS.

This establishes an empirical within-run variability scale without introducing a detector model.

## Day-to-day versus within-run variation

At each anchor the diagnostic compares the absolute day-to-day ASD ratio in dB against the maximum absolute block/full-repeat ASD ratio in dB.

It also compares the full 1--200 kHz RMS day-to-day shape difference against the largest block RMS.

Important output fields are:

```text
day_vs_within_repeat.n_anchors_day_exceeds_all_blocks
day_vs_within_repeat.fraction_anchors_day_exceeds_all_blocks
day_vs_within_repeat.day_rms_over_block_rms_max_ratio
day_vs_within_repeat.day_fit_band_rms_exceeds_all_blocks
```

If the day-to-day shape exceeds every same-run block fluctuation, the spectral difference is unlikely to be explained solely by finite-record statistics within the repeat run.

This is not a formal stationarity confidence interval. The blocks are deterministic contiguous subsets, not independent experimental repetitions.

## Interpretation flags

### `day_shape_exceeds_within_repeat_block_variation`

True when the 1--200 kHz RMS day-to-day shape difference exceeds the largest within-repeat block RMS.

### `transfer_ratio_tracks_day_shape_broadly`

A deliberately simple diagnostic screen requiring:

- dB-shape correlation >= 0.8;
- empirical-minus-transfer RMS <= 1 dB across 1--200 kHz.

This is a screening flag, not a statistical hypothesis test.

### `strong_electronics_drift_claim_allowed`

Always false.

Even if the empirical ASD ratio closely follows the fitted transfer ratio, the repeat detector operating point was not independently calibrated/refit. The result therefore supports a measurement-chain-like drift interpretation but does not uniquely identify electronics drift.

## Run

From the repository root:

```powershell
python PoST_Simulations/subScript/readout_day_to_day_drift_diagnostic.py
```

Default output:

```text
PoST_Simulations/.noise_optimization_work_rsh_sweep/readout_day_to_day_drift_diagnostic.json
```

## What to look at next

1. **Day-to-day ratio is larger than all block variation, and the transfer ratio matches it well.** This strengthens the interpretation that the fitted transfer drift corresponds to a reproducible measurement-chain-like spectral change.
2. **Day-to-day ratio is larger than block variation, but the transfer ratio does not match it.** The local readout fit is likely absorbing other run-dependent physics or operating-point differences.
3. **Day-to-day ratio is comparable to block variation.** The apparent transfer drift is not distinguishable from within-run variability at the current record count.
4. **Only narrow frequencies exceed block variation.** Inspect line/comb contamination before assigning a broadband transfer interpretation.
