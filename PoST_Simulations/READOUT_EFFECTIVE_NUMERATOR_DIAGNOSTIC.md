# Effective numerator / minimal-order diagnostic

## Purpose

The low/mid identifiability profiles showed that the 1--40 kHz data constrain a pole-like scale near 12--13 kHz much more robustly than they constrain individual zero and lead-zero corner frequencies.

Rather than continuing to interpret degenerate corner locations, this diagnostic asks:

> How many independent numerator-shape degrees of freedom are actually required by the 1--40 kHz data?

## Nested model family

All models use the same second-order pole denominator.

Let:

```text
x = f / 40 kHz
P(x) = numerator magnitude squared
```

The fitted complexity ladder is strictly nested:

### Order 0: pole only

```text
P(x) = 1
```

Total free parameters: 2

- pole frequency
- pole Q

### Order 1: one numerator degree

```text
P(x) = 1 + v^2 x^2
```

Total free parameters: 3.

`v=0` recovers order 0 exactly.

### Order 2: two numerator degrees

```text
P(x) = (1 + u x^2)^2 + v^2 x^2
```

Total free parameters: 4.

`u=0` recovers order 1 exactly.

### Order 3: three numerator degrees

```text
P(x) = [(1 + u x^2)^2 + v^2 x^2] (1 + w^2 x^2)
```

Total free parameters: 5.

`w=0` recovers order 2 exactly.

Order 3 is algebraically equivalent to the numerator of the current five-parameter infinity-lead-pole hybrid transfer. It is therefore a reparameterization, not a new more-flexible model.

## Polynomial coefficients

Every fitted numerator is reported as:

```text
P(x) = 1 + c2 x^2 + c4 x^4 + c6 x^6
```

The output reports the dimensionless coefficients `c2`, `c4`, and `c6`.

These coefficients describe effective magnitude shape only. They are not resistor, capacitor, pole, or zero component values.

## Why use the factorized form?

A directly fitted polynomial can become negative at frequencies outside the 1--40 kHz fit band, which is invalid for a magnitude-squared transfer and can create problems in the first-alias branch.

The nested factorization above guarantees:

```text
P(x) > 0
```

for all frequencies while still providing a 0/1/2/3 numerator-degree ladder.

## Order-3 equivalence check

The tracked low/mid free hybrid fit is analytically mapped into the order-3 effective parameterization.

At runtime the diagnostic compares the normalized model produced by the mapped order-3 coefficients against the original hybrid model.

The run aborts if:

```text
max absolute normalized-model difference > 1e-10
```

This ensures that any difference between complexity orders comes from removing degrees of freedom, not from changing the five-parameter family.

## Git-tracked inputs

Configuration:

```text
PoST_Simulations/config/readout_effective_numerator_config.json
```

Tracked low/mid free hybrid snapshot:

```text
PoST_Simulations/cases/tagawa_20241205_r1ch12_215mK_1400uA_gain5_repeat/lowmid_free_transfer_snapshot.json
```

Shared reference transfer:

```text
PoST_Simulations/config/shared_readout_reference_transfer.json
```

The post-filter white floor remains fixed to the preceding profiled value.

## Fit and holdout

Optimization uses only:

```text
1 kHz <= f < 40 kHz
```

The 40--200 kHz region is evaluated only after each order is fitted:

```text
40 kHz <= f <= 200 kHz
```

No holdout point is passed to the optimizer.

## Optimizer candidates

For each order the diagnostic considers:

- differential evolution;
- DE -> least-squares refinement;
- exact shared-transfer warm start mapped into that order;
- shared warm start -> least-squares;
- exact preceding free-fit warm start mapped into that order;
- free-fit warm start -> least-squares.

The candidate with the lowest actual production shape score is retained.

## Minimal-order screen

An order is considered `near_free_hybrid` only if both are satisfied relative to the tracked five-parameter free hybrid:

```text
score / free_hybrid_score <= 5
RMS residual increase <= 0.05 dB
```

The lowest order passing both conditions is reported as:

```text
minimal_order_near_free_hybrid
numerator_dof_required_to_match_free
```

Important flags are:

```text
interpretation_flags.pole_only_sufficient
interpretation_flags.one_numerator_degree_sufficient
interpretation_flags.two_numerator_degrees_sufficient
interpretation_flags.three_numerator_degrees_required
```

## Incremental complexity gain

For each step 0->1, 1->2, and 2->3 the output reports:

- new/old fit score ratio;
- RMS improvement in dB;
- new/old high-frequency holdout score ratio.

This helps distinguish a genuinely necessary extra shape degree from a tiny in-band improvement that only increases flexibility.

## Interpretation

If order 0 is sufficient, the low/mid data require only a pole section.

If order 1 is sufficient, the data require one effective numerator curvature degree in addition to the pole.

If order 2 is sufficient but order 1 is not, two numerator coefficients are needed, but the third hybrid numerator degree is not identified by the low/mid data.

If only order 3 reaches the near-free screen, all three numerator shape degrees are required to reproduce the five-parameter hybrid fit over 1--40 kHz.

Even in that last case, the coefficients should be interpreted as effective magnitude-shape degrees rather than individual physical zero/lead-zero corners.

## Run

From the repository root:

```powershell
python PoST_Simulations/subScript/readout_effective_numerator_diagnostic.py
```

Default output:

```text
PoST_Simulations/.noise_optimization_work_rsh_sweep/readout_effective_numerator_diagnostic.json
```
