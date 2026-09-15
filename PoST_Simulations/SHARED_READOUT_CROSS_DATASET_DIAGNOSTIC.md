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

## Git-tracked inputs

The standard diagnostic no longer requires user-created JSON files outside the
repository.  Its required configuration is tracked in Git.

Default manifest:

```text
PoST_Simulations/config/shared_readout_cross_dataset_manifest.json
```

Default shared-transfer reference:

```text
PoST_Simulations/config/shared_readout_reference_transfer.json
```

The reference detector/fit snapshot currently used by the manifest is also
tracked with the case:

```text
PoST_Simulations/cases/tagawa_20241206_r1ch12_215mK_1400uA_gain5_day2/shared_readout_detector_snapshot.json
```

That compact snapshot intentionally stores only the fields required by this
diagnostic: the fit definition and frozen `best_case_parameters`.  This avoids
depending on a separate local optimizer-summary file whose contents may drift
from the diagnostic result being cross-validated.

The existing case `comparison_summary.json` is already Git tracked and is
referenced directly by the manifest.

When a new same-readout validation dataset is added, add its detector snapshot
and comparison summary under `PoST_Simulations/cases/`, then add a
`role: "validation"` entry to the tracked manifest.  The raw CH0 records may
remain on the acquisition drive; `experiment_path` can stay in the tracked
comparison summary or be overridden in the manifest if the mount point differs.

Relative `summary` and `comparison_summary` paths are resolved relative to
the tracked manifest location.

### Roles

Use:

- `reference` for the dataset from which the shared transfer was originally
  inferred;
- `validation` for independent bias, bath, or detector conditions that use
  the same readout chain.

A manifest containing only reference cases can exercise the code, but the
diagnostic explicitly refuses to call that cross-validation.

## Reference transfer

The current accepted transfer is frozen in the tracked file:

```text
PoST_Simulations/config/shared_readout_reference_transfer.json
```

It records the exact best infinity-pole result from the latest pole-profile
diagnostic, including provenance and the source shape score.  The diagnostic
reads its `best_profile_row` and therefore preserves the exact infinite
lead/lag-pole state rather than approximating it with a multi-MHz finite pole.

A later diagnostic result should update this tracked reference file in the same
commit series that updates the cross-dataset assumptions.

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

With the tracked defaults, run from the repository root with no input-JSON
arguments:

```powershell
python PoST_Simulations/subScript/shared_readout_cross_dataset_diagnostic.py
```

Default output:

```text
PoST_Simulations/.noise_optimization_work_rsh_sweep/shared_readout_cross_dataset_diagnostic.json
```

For a cheaper fixed-transfer-only run:

```powershell
python PoST_Simulations/subScript/shared_readout_cross_dataset_diagnostic.py --skip-local-refit
```

`--manifest` and `--reference-profile` remain available only as explicit
overrides for diagnostic experiments.  They are not required for the standard
Git-reproducible workflow.

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
