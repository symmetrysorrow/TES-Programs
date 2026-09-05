# Conditional proxy ensemble comparison

Strict target conclusion: **C — exact target physical case remains unidentified**.

Exploratory conclusion: **P4 — proxy parameter space is still too underconstrained to make a useful target reproduction statement**.

The Stage-A range was frozen before this comparison. Sampled q05/q95 are descriptive quantiles, not probabilities.

## Pre-analysis shape

| Hz | experiment | proxy min | proxy q05 | proxy q50 | proxy q95 | proxy max | status |
|---:|---:|---:|---:|---:|---:|---:|---|
| 10 | 46.1562 | 0.999258 | 1.02361 | 1.23326 | 1.70886 | 2.12984 | outside_sampled_envelope |
| 100 | 4.02433 | 1.00608 | 1.02978 | 1.23352 | 1.68924 | 2.05604 | outside_sampled_envelope |
| 1000 | 1 | 1 | 1 | 1 | 1 | 1 | inside_sampled_min_max |
| 3000 | 0.935658 | 0.541824 | 0.590047 | 0.696302 | 0.840284 | 0.954697 | inside_sampled_min_max |
| 5000 | 0.862552 | 0.430532 | 0.477168 | 0.598108 | 0.80256 | 0.882456 | inside_sampled_min_max |
| 7000 | 0.782724 | 0.391046 | 0.421822 | 0.556263 | 0.790852 | 0.806565 | inside_sampled_min_max |
| 10000 | 0.636483 | 0.348106 | 0.389138 | 0.515133 | 0.785704 | 0.803826 | inside_sampled_min_max |

## Band diagnostics

- **10-100_Hz**: partially_or_fully_outside (2 sampled-anchor points outside min/max).
- **100-1000_Hz**: partially_or_fully_outside (1 sampled-anchor points outside min/max).
- **1-3_kHz**: covered_by_sampled_min_max (0 sampled-anchor points outside min/max).
- **3-10_kHz**: covered_by_sampled_min_max (0 sampled-anchor points outside min/max).

Post-analysis comparison applies the same deterministic 2nd-order 10 kHz Bessel zero-phase ASD factor. It is secondary; detector-side shape interpretation uses pre-analysis.

No parameter range was adjusted and no best-looking member was promoted to a parameter estimate.
