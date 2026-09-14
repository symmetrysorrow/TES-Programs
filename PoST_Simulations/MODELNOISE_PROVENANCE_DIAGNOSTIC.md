# modelnoise.txt provenance diagnostic

## Purpose

The empirical digital-transfer diagnostic showed that the measured 10 kHz
analysis transfer agrees with the production second-order Bessel
`filtfilt` response to far below the scale of the detector residual. It also
showed that the stored `CH0_noise/modelnoise.txt` shape is not identical to a
fresh canonical reconstruction from the same raw data.

This diagnostic isolates that target-file provenance question before any more
detector-physics fitting.

It never overwrites `modelnoise.txt`, changes no production code, and keeps
the detector parameters fixed.

## Inputs

The script requires:

- the target-case `comparison_summary.json`, which supplies the exact
  accepted CH0 record indices and acquisition settings;
- the optimizer `summary.json`, which supplies the frozen best detector
  parameters and fit-objective settings.

The experiment root is normally read from `comparison_summary.json` and can
be overridden with `--experiment-path`.

## Fresh canonical target

The fresh target is reconstructed from exactly the accepted record mask already
stored in `comparison_summary.json`. No record-selection decision is rerun.

Each accepted raw record follows the current canonical pipeline:

```text
mean removal
-> 10 kHz second-order digital Bessel filtfilt
-> Hann window
-> rFFT power
-> power average across records
-> one-sided ASD
```

The absolute voltage/current calibration is irrelevant to the primary
comparison because stored and fresh targets are independently normalized at
1 kHz.

## Stored/fresh shape audit

The diagnostic computes

```text
20 log10(stored_normalized / fresh_normalized)
```

from 1 to 200 kHz.

To distinguish broad shape changes from narrow spectral features, the dB
difference is split into:

- a broad component using a 1 kHz median trend by default;
- a narrow component equal to total minus that trend.

The JSON records total, broad, and narrow RMS/max values, band summaries, the
fraction of bins with narrow excursions above 0.5 dB, and up to 12 separated
largest narrow excursions.

The smoothing width is a diagnostic scale, not a physical filter assumption,
and can be changed with `--smooth-width-hz`.

## Frozen detector comparison

The same frozen best detector spectrum is evaluated against two targets on the
same optimizer fit grid:

1. the stored `modelnoise.txt`, using the historical frequency mapping used
   by `Opt_noise.target_spectrum`;
2. the fresh canonical raw-record reconstruction.

Both use the same production `fit_score` and band diagnostics.

The output reports:

- stored and fresh shape scores;
- whether the stored score reproduces `summary.json["best_case_score"]`;
- the score change caused only by replacing the target file;
- 1--5, 5--15, 15--40, 40--100, and 100--200 kHz model/target residuals;
- whether the characteristic 5--15 kHz deficit plus 40--100 kHz excess remains
  for the fresh target.

A target replacement is performed only in memory for the diagnostic. The
stored target file is never modified.

## Run

From the repository root in PowerShell:

```powershell
python PoST_Simulations/subScript/modelnoise_provenance_diagnostic.py --comparison-summary PoST_Simulations/cases/tagawa_20241206_r1ch12_215mK_1400uA_gain5_day2/comparison_summary.json --optimizer-summary PoST_Simulations/.noise_optimization_work_rsh_sweep/summary.json
```

Default output:

```text
PoST_Simulations/.noise_optimization_work_rsh_sweep/modelnoise_provenance_diagnostic.json
```

If the raw-data drive is mounted elsewhere:

```powershell
python PoST_Simulations/subScript/modelnoise_provenance_diagnostic.py --comparison-summary PoST_Simulations/cases/tagawa_20241206_r1ch12_215mK_1400uA_gain5_day2/comparison_summary.json --optimizer-summary PoST_Simulations/.noise_optimization_work_rsh_sweep/summary.json --experiment-path G:/tagawa/20241206/r1ch12_215mK_1400uA1400uA_difftrig5e-5_rate500k_samples100k_gain5_day2
```

## Interpretation

If the fresh target retains the 5--15 kHz model deficit and 40--100 kHz model
excess with similar magnitude, the stored-target provenance mismatch is not the
main cause of the detector residual. The investigation can then move to the
analog/readout hardware response with the digital analysis stage independently
validated.

If the residual pattern changes strongly or disappears for the fresh target,
the target-file generation history must be resolved before further detector
parameter interpretation.

A large narrow-component maximum together with a small broad-component RMS
would indicate isolated lines or spikes rather than a broad transfer-shape
problem. A substantial broad component would instead show that stored and
fresh targets differ systematically across frequency.
