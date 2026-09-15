# Hybrid readout lead/lag pole profile and alias-identifiability diagnostic

## Purpose

The matched-topology comparison found that the best six-parameter pre-ADC
hybrid combines:

- a general second-order pole/zero section;
- one real lead/lag pole/zero pair.

Its lead/lag pole landed at approximately the previous 300 kHz search ceiling.
That boundary hit is ambiguous. It can mean either:

1. a real finite high-frequency hardware corner lies near or above the search
   limit; or
2. the fit only wants the limiting shape obtained as that pole moves toward
   infinity, so the pole itself is not identifiable.

This diagnostic profiles that one parameter explicitly before adding any more
model complexity.

## Frozen physics and target

The detector model remains frozen. The target is the fresh pre-analysis ASD
reconstructed from the exact accepted CH0 record indices in
`comparison_summary.json`.

The model path remains:

```text
frozen intrinsic detector ASD
-> production 100 kHz order-4 mag-normalized analog Bessel
-> diagnostic hybrid pre-ADC transfer
-> first ADC alias fold
-> production post-filter white term
-> normalized pre-analysis ASD
```

No 10 kHz digital analysis filter is applied.

## Profile definition

The lead/lag pole is fixed at each default point:

```text
150 kHz
250 kHz
300 kHz
500 kHz
1 MHz
2 MHz
5 MHz
infinity
```

At every fixed pole, exactly five parameters are refit:

- second-order pole center;
- second-order pole Q;
- second-order zero center;
- second-order zero Q;
- lead/lag zero corner.

The second-order Q range remains 0.10--20, so overdamped and complex-conjugate
solutions are both allowed.

The exact infinity point is not approximated by a very large number. It removes
the lead/lag denominator analytically:

```text
sqrt(1 + (f / f_zero)^2)
```

This is the mathematical `f_pole -> infinity` limit.

## Why the infinity point matters

If the score keeps improving as the fixed pole is moved upward and the exact
infinity row is equal to or better than the 5 MHz row, the data do not identify
a finite high-frequency pole. In that case, interpreting a fitted value such as
300 kHz or several MHz as a physical hardware corner would be misleading.

The useful information would instead be that the measured shape requires a
broad rising lead/shelf contribution over the fit band.

If the score has a clear finite minimum and worsens substantially again before
the infinity row, a finite high-frequency corner becomes a more plausible
readout feature.

The JSON reports both conditions through:

- `supports_pole_to_infinity_nonidentifiability`;
- `supports_finite_high_frequency_corner`.

These are diagnostic flags, not statistical confidence intervals.

## Alias bookkeeping

For the production baseline and for the best profile row, the diagnostic also
reports PSD fractions at:

```text
1, 5, 10, 20, 40, 70, 100, 150, 200 kHz
```

Each row contains:

- main-branch PSD fraction;
- first-alias PSD fraction;
- post-filter-white PSD fraction;
- first-alias / main ASD ratio.

The fractions are evaluated before the final 1 kHz normalization, so they are
true bookkeeping fractions within this first-alias model.

This makes it possible to determine whether the high-pole profile is changing
the fit mainly by modifying the measured-frequency branch or by changing the
folded contribution from frequencies near `rate - f`.

## Optimization

Each fixed-pole row uses:

1. differential evolution;
2. bounded least-squares polish.

The production shape objective and broad-band penalties are reused unchanged.

Default center/corner search range:

```text
1 kHz -- 300 kHz
```

Default second-order Q search range:

```text
0.10 -- 20
```

The fixed lead/lag pole itself is allowed above Nyquist because it is an analog
pre-ADC transfer parameter. The transfer is evaluated at both the main
frequency and first-alias frequency before folding.

## Run

From the repository root:

```powershell
python PoST_Simulations/subScript/preanalysis_hybrid_pole_profile_diagnostic.py --summary PoST_Simulations/.noise_optimization_work_rsh_sweep/summary.json --comparison-summary PoST_Simulations/cases/tagawa_20241206_r1ch12_215mK_1400uA_gain5_day2/comparison_summary.json
```

Default output:

```text
PoST_Simulations/.noise_optimization_work_rsh_sweep/preanalysis_hybrid_pole_profile_diagnostic.json
```

A custom finite profile can be supplied with:

```powershell
--pole-grid-hz 150000,250000,300000,500000,1000000,2000000,5000000
```

The exact infinity row is always added automatically.

## Important output fields

### `profile_rows`

One refit result for every fixed finite pole plus the exact infinity limit.

Each row records:

- shape score;
- score ratio to the frozen baseline;
- refitted parameters;
- root regimes and equivalent real roots when Q <= 0.5;
- broad-band residuals;
- numerical boundary positions;
- the 1 dB RMS / 3 dB maximum-residual screen.

### `profile_trend`

Summarizes the identifiability question:

- best finite pole;
- best finite score;
- exact infinity score;
- infinity / best-finite score ratio;
- whether the score is non-increasing from 300 kHz to 5 MHz;
- whether infinity is the best profile point;
- finite-corner and pole-to-infinity interpretation flags.

### Alias fractions

`baseline.alias_fraction_sample` and
`best_profile_alias_fraction_sample` make the role of first-alias power
explicit.

## Interpretation

### Pole-to-infinity supported

If the profile improves monotonically above 300 kHz and the infinity row is
within roughly one percent of or better than the 5 MHz row, the high-frequency
pole is not identified.

Do not translate the best finite value into a component time constant. The
meaningful fitted feature is the lead/lag zero and the resulting rising shelf.

The next step should then be hardware/schematic comparison for a gain slope,
zero, differentiating stage, or a more complete anti-alias/readout transfer,
rather than searching for a multi-MHz pole.

### Finite corner supported

If a finite pole well below the top of the grid gives a distinctly better score
than both higher poles and the infinity limit, that corner is worth comparing
against readout electronics, anti-alias stages, cabling, amplifier bandwidth,
or SQUID/input-network constants.

### Alias fraction changes strongly

If the preferred profile changes the fit largely through the first-alias PSD
fraction, the hardware transfer above Nyquist is materially coupled to the
in-band result. A single first-alias approximation should then be checked
against a wider alias sum before assigning detailed circuit meaning.

### Alias fraction remains small

If first-alias power remains negligible, the fitted profile is primarily an
in-band transfer-shape result. That would make the inferred ~5 kHz zero and
~10--100 kHz second-order structure more direct targets for hardware
comparison.

## Guardrail

A pole that runs toward infinity is a non-identifying mathematical limit, not
evidence for a physical multi-MHz component. Likewise, first-alias fractions
are bookkeeping within the current sampling model; they do not substitute for a
directly measured analog transfer function.
