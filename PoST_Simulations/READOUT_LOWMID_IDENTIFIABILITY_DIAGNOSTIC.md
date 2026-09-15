# Low/mid transfer identifiability profiles

## Purpose

The strict 1--40 kHz holdout fit produced an excellent low/mid match but drove several hybrid-transfer parameters to very different values, including an overdamped zero section with equivalent roots near 4.36 kHz and 204 kHz plus a lead zero near 81 kHz.

That result is a warning that the five-parameter transfer may be non-identifiable over the restricted fit band.

This diagnostic therefore asks a different question:

> Which transfer parameters can actually be constrained by the 1--40 kHz data after the white floor is fixed?

## Git-tracked inputs

Configuration:

```text
PoST_Simulations/config/readout_lowmid_identifiability_config.json
```

Tracked low/mid free-fit snapshot:

```text
PoST_Simulations/cases/tagawa_20241205_r1ch12_215mK_1400uA_gain5_repeat/lowmid_free_transfer_snapshot.json
```

Reference shared transfer:

```text
PoST_Simulations/config/shared_readout_reference_transfer.json
```

The low/mid snapshot freezes the preceding strict-holdout free optimum and the profiled white floor. It is a numerical reference only, not a hardware calibration.

## Fixed-parameter profiles

For each profile point one selected transfer parameter is held fixed while all remaining hybrid parameters are re-optimized.

The default profiles are:

```text
zero_Hz:
30, 40, 50, 60, 68.5(shared), 80, 100, 150, 300 kHz

leadlag_zero_Hz:
3, 4, 4.94(shared), 6, 8, 10, 20, 50, 80 kHz

pole_Hz:
9, 10, 11, 11.82(free), 13, 13.34(shared), 15, 18 kHz
```

Each profile point uses:

- differential evolution;
- DE -> least-squares refinement;
- exact shared warm start;
- shared -> least-squares warm start;
- exact preceding free-fit warm start;
- free-fit -> least-squares warm start.

The lowest actual score among these candidates is retained.

## Near-free screen

A fixed profile point is marked `near_free` only when both are true:

```text
score / free_score <= 5
RMS residual increase <= 0.05 dB
```

The score-ratio threshold is intentionally paired with an RMS threshold because the free-fit score is extremely small; a modest absolute degradation can otherwise look artificially large as a ratio.

A profile is flagged as broadly non-identified when at least three near-free points span a factor of at least two in the fixed parameter.

Important fields:

```text
fixed_parameter_profiles
profile_summaries
profile_summaries.<parameter>.near_free_span_ratio
profile_summaries.<parameter>.broad_nonidentifiability_over_tested_range
```

## Reduced-topology fits

The diagnostic also asks whether the extreme free-fit values are necessary at all.

Three reduced fits are run:

### `shared_zero_section_fixed`

Fix:

```text
zero_Hz = shared value
zero_Q  = shared value
```

Refit only pole frequency, pole Q, and lead zero.

### `shared_lead_zero_fixed`

Fix the lead zero to the shared value and refit the pole and second-order zero section.

### `shared_zero_section_and_lead_zero_fixed`

Fix the shared zero section and shared lead zero, leaving only pole frequency and pole Q free.

Each reduced fit is compared with the five-parameter free low/mid optimum using the same near-free screen.

Important fields:

```text
reduced_topology_fits
interpretation_flags.shared_zero_section_reduced_fit_near_free
interpretation_flags.shared_lead_zero_reduced_fit_near_free
interpretation_flags.shared_zero_and_lead_reduced_fit_near_free
```

## Interpretation

If `zero_Hz` remains near-free across a broad range, neither ~75 kHz nor ~30 kHz should be interpreted as an identified zero frequency.

If `leadlag_zero_Hz` is similarly flat, the 81 kHz lead zero from the restricted free fit is also only a compensating parameter.

If the pole-frequency profile is much sharper than the zero/lead profiles, that would support the interpretation that the robust low/mid information is concentrated near the ~10--15 kHz pole scale while the other sections are weakly identified.

If the shared zero section can be fixed with little RMS penalty, the extreme overdamped zero from the free low/mid fit is unnecessary.

If even the two-parameter reduced fit with shared zero and lead sections remains near-free, the low/mid data are effectively constraining only the pole section.

Conversely, a sharp profile or large reduced-fit penalty identifies which additional degree of freedom the low/mid data genuinely require.

None of these outcomes identifies a physical electronics component by itself.

## Run

From the repository root:

```powershell
python PoST_Simulations/subScript/readout_lowmid_identifiability_diagnostic.py
```

Default output:

```text
PoST_Simulations/.noise_optimization_work_rsh_sweep/readout_lowmid_identifiability_diagnostic.json
```
