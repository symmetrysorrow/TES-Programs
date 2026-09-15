# Shared readout cross-dataset validation diagnostic

## Purpose

The current pre-analysis diagnostics identify a compact readout-transfer shape
that strongly improves the reference dataset while detector physics is frozen:

- a real lead zero near 5 kHz;
- a second-order pole near 13 kHz;
- a second-order zero near 68 kHz;
- the lead/lag pole running to the exact infinity limit.

That result is still derived from one detector dataset.  The decisive next test
is whether the same transfer works on independent detector conditions acquired
through the same readout chain.

This diagnostic performs that cross-dataset test without refitting the shared
readout transfer.

## What is fixed and what changes

Across every case:

### Held fixed globally

The shared pre-ADC readout transfer is taken from the
`best_profile_row` of a
`preanalysis_hybrid_pole_profile_diagnostic.json` file.

The transfer parameters are not refit for individual cases.

For the current reference result this corresponds to the infinity-pole hybrid:

```text
general second-order pole / zero
x
real lead zero with lead pole -> infinity
```

### Case specific but frozen

Each manifest case supplies its own optimizer `summary.json`.

The diagnostic reads that case's `best_case_parameters` and freezes those
detector/source parameters for that case.

Thus the cross-dataset comparison changes detector conditions from case to
case, but does not allow the shared readout transfer to move.

### Experimental target

Each case also supplies a `comparison_summary.json`.

The experimental target is freshly reconstructed from that case's exact
`accepted_record_indices`:

```text
raw CH0 records
-> exact stored accepted mask
-> mean removal
-> Hann
-> power average
-> one-sided ASD
```

No 10 kHz digital analysis filter is applied to target or model.

## Model path

For every case:

```text
case-specific frozen intrinsic detector ASD
-> production 100 kHz order-4 mag-normalized analog Bessel
-> one globally shared diagnostic pre-ADC transfer
-> first alias fold
-> case-specific production post-filter white term
-> 1 kHz normalization
```

The diagnostic checks that the acquisition rate in the case comparison summary
matches the rate used by the case-specific detector summary.

## Manifest

The diagnostic accepts one JSON manifest so local datasets that are not
committed to the repository can be included directly.

Example:

```json
{
  "cases": [
    {
      "label": "215mK_1400uA_reference",
      "role": "reference",
      "summary": "D:/work/215mK_1400uA/summary.json",
      "comparison_summary": "D:/Github/TES-Programs/PoST_Simulations/cases/tagawa_20241206_r1ch12_215mK_1400uA_gain5_day2/comparison_summary.json"
    },
    {
      "label": "same_readout_other_bias",
      "role": "validation",
      "summary": "D:/work/other_bias/summary.json",
      "comparison_summary": "D:/work/other_bias/comparison_summary.json",
      "experiment_path": "G:/tagawa/20241206/<raw-data-directory>"
    }
  ]
}
```

`experiment_path` is optional.  When omitted, it is read from that case's
`comparison_summary.json`.

Relative `summary` and `comparison_summary` paths are resolved relative to
the manifest location.

### Roles

Use:

- `reference` for the dataset from which the shared transfer was originally
  inferred;
- `validation` for independent bias, bath, or detector conditions that use
  the same readout chain.

A manifest containing only reference cases can exercise the code, but the
diagnostic explicitly refuses to call that cross-validation.

## Reference profile

Pass the latest pole-profile JSON:

```text
PoST_Simulations/.noise_optimization_work_rsh_sweep/preanalysis_hybrid_pole_profile_diagnostic.json
```

The diagnostic reads `best_profile_row` and preserves the exact finite or
infinite lead/lag-pole state.

For the present result, the infinity limit is therefore reused exactly rather
than approximated with a multi-MHz finite pole.

## Shared fixed evaluation

For every case the output contains:

- production baseline shape score;
- fixed-shared-transfer shape score;
- shared / baseline score ratio;
- broad-band residuals;
- change in absolute broad-band mean error;
- RMS / maximum dB residuals;
- first-alias PSD fractions;
- representative 1--200 kHz curves.

The main broad-band flags are:

- `shared_improves_shape_score`;
- `shared_improves_5_15k`;
- `shared_improves_40_100k`;
- `shared_does_not_worsen_100_200k_mean`;
- `shared_strict_broad_improvement`.

## Local-refit comparator

By default the diagnostic also performs a separate local fit for each case
using the same topology and the same fixed lead/lag-pole state.

This local fit is **not** used as the shared cross-validation model.

Its only purpose is to ask:

> How much fit quality is lost by forcing the reference transfer to remain
> shared?

The output reports:

- local best score;
- shared / local-best score ratio;
- local parameter drift relative to the shared transfer;
- RMS and maximum log10 parameter-ratio drift.

The default diagnostic tolerance is:

```text
shared score <= 1.25 x local best score
```

This is a pragmatic diagnostic tolerance, not a statistical confidence level.

Use `--skip-local-refit` when only the strict fixed-transfer evaluation is
wanted.

## Aggregate interpretation

The aggregate output separately counts reference and validation cases.

The flag:

```text
supports_shared_readout_transfer_across_validation_cases
```

requires all validation cases to:

1. improve total shape score with the shared transfer;
2. improve both the 5--15 kHz and 40--100 kHz broad residuals;
3. when local refits are enabled, keep the shared score within the configured
   shared/local tolerance.

If no validation case is present:

```text
validation_dataset_required_for_cross_validation_claim = true
supports_shared_readout_transfer_across_validation_cases = false
```

This prevents the reference dataset from validating itself.

## Run

Create a manifest, for example:

```text
D:/work/shared_readout_manifest.json
```

Then run from the repository root:

```powershell
python PoST_Simulations/subScript/shared_readout_cross_dataset_diagnostic.py --manifest D:/work/shared_readout_manifest.json --reference-profile PoST_Simulations/.noise_optimization_work_rsh_sweep/preanalysis_hybrid_pole_profile_diagnostic.json
```

Default output is written beside the manifest:

```text
shared_readout_cross_dataset_diagnostic.json
```

For a cheaper fixed-transfer-only run:

```powershell
python PoST_Simulations/subScript/shared_readout_cross_dataset_diagnostic.py --manifest D:/work/shared_readout_manifest.json --reference-profile PoST_Simulations/.noise_optimization_work_rsh_sweep/preanalysis_hybrid_pole_profile_diagnostic.json --skip-local-refit
```

## How to interpret the result

### Strong support for a shared readout transfer

The readout interpretation becomes substantially stronger when independent
validation cases show all of the following:

- the fixed reference transfer improves their baseline scores;
- the characteristic 5--15 kHz deficit and 40--100 kHz excess both shrink;
- the shared transfer performs close to each case's local best transfer;
- local transfer parameters do not drift strongly from the reference values.

In that situation, changing detector operating conditions while preserving the
same readout transfer is exactly the behavior expected for a measurement-chain
effect.

### Fixed shared transfer improves, but local parameters drift

This is mixed evidence.

The generic transfer topology has useful shape leverage across datasets, but
the actual frequency scales may be absorbing detector-dependent electrical
dynamics.  Production code should not yet treat the transfer as a fixed
electronics response.

### Each local fit works but the shared fixed transfer fails

That result argues against a single readout-chain transfer as the explanation.
The fitted transfer would be functioning mainly as a flexible proxy for
case-dependent missing physics.

### Reference case only

This is not cross-validation, regardless of how good the fit is.

## Guardrail

This diagnostic tests consistency of one magnitude-only pre-ADC transfer
across datasets.  Even a successful cross-dataset result does not by itself
identify the physical circuit or its phase.  Hardware schematics, direct
transfer measurements, or controlled electronics calibration remain necessary
before moving the diagnostic transfer into the production noise model.
