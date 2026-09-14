# Empirical digital-transfer diagnostic

## Purpose

The measurement-chain convention sweep showed that changing only known
measurement-chain conventions can improve the mid/high-band residuals while
all detector parameters remain frozen. The strongest numerical proxy used a
10 kHz fourth-order single-pass digital Bessel response, but the production
experimental code explicitly applies a second-order digital Bessel with
`scipy.signal.filtfilt`.

This diagnostic therefore measures the digital transfer directly from the
experimental raw CH0 records instead of inferring it from the TES noise fit.

No TES model is called.

## Inputs and record selection

The script takes the target-case `comparison_summary.json`. It reuses the
exact `accepted_record_indices` stored there and does not rerun or alter the
record-selection logic.

For those same accepted raw records it computes:

1. pre-analysis ASD: mean removal -> Hann -> rFFT -> power average;
2. post-analysis ASD: mean removal -> 10 kHz second-order Bessel
   `filtfilt` -> Hann -> rFFT -> power average.

The empirical transfer is

```text
H_empirical(f) = ASD_post(f) / ASD_pre(f)
```

Because both spectra use the same raw voltage records, this ratio is
dimensionless and does not require the unresolved voltage-to-current
calibration.

## Transfer comparisons

The empirical transfer is compared with three predictions:

- production analytic response: second-order digital Bessel, two magnitude
  passes, `|H_2|^2`;
- fourth-order single-pass proxy, `|H_4|`, included only because that shape
  was favored by the preceding convention screen;
- paired finite-record production simulation.

The paired finite-record simulation synthesizes Gaussian records using the
measured pre-analysis ASD as the input spectrum. Each synthetic record is
estimated both before and after the production filter, with the same Hann
window and one-sided ASD estimator. This preserves finite-record, edge,
window, and spectral-leakage effects while avoiding TES physics.

The JSON reports absolute and independently 1 kHz-normalized residuals in dB,
including the same 1--5, 5--15, 15--40, 40--100, and 100--200 kHz bands used
throughout the noise investigation.

## Stored modelnoise audit

The script also loads

```text
<experiment>/CH0_noise/modelnoise.txt
```

when present and compares its normalized shape with the newly recomputed
canonical post-analysis ASD.

This check is intentionally insensitive to the absolute eta calibration. A
shape mismatch would indicate that the optimizer target file was generated
under different processing/data-selection semantics than the current
canonical raw-record pipeline.

## Run

From the repository root in PowerShell:

```powershell
python PoST_Simulations/subScript/empirical_digital_transfer_diagnostic.py --comparison-summary PoST_Simulations/cases/tagawa_20241206_r1ch12_215mK_1400uA_gain5_day2/comparison_summary.json
```

The experiment path is normally read from `comparison_summary.json`. It can
be overridden when the raw-data drive is mounted elsewhere:

```powershell
python PoST_Simulations/subScript/empirical_digital_transfer_diagnostic.py --comparison-summary PoST_Simulations/cases/tagawa_20241206_r1ch12_215mK_1400uA_gain5_day2/comparison_summary.json --experiment-path G:/tagawa/20241206/r1ch12_215mK_1400uA1400uA_difftrig5e-5_rate500k_samples100k_gain5_day2
```

Default output:

```text
PoST_Simulations/cases/tagawa_20241206_r1ch12_215mK_1400uA_gain5_day2/empirical_digital_transfer_diagnostic.json
```

By default the paired finite-record simulation uses the same number of records
as the accepted experimental set. `--finite-records` can override that
count, and `--finite-seed` controls its deterministic RNG seed.

## Interpretation

If the empirical post/pre transfer is best reproduced by the production
finite-record prediction, the current second-order `filtfilt` convention is
supported and the fourth-order/single-pass numerical improvement should be
treated as a proxy for some other missing shape.

If the empirical transfer is closer to the fourth-order/single-pass proxy, or
if the stored `modelnoise.txt` shape disagrees with a fresh canonical
reconstruction, the experimental processing/target-file provenance must be
resolved before further detector-physics fitting.

This diagnostic does not test the 100 kHz analog hardware response and cannot
identify TES physics.
