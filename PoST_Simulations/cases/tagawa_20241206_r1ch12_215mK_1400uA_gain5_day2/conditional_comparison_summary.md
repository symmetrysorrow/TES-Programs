# Conditional proxy ensemble comparison

Strict target conclusion: **C — exact target physical case remains unidentified**.

Exploratory conclusion: **P4 — proxy parameter space is still too underconstrained to make a useful target reproduction statement**.

The Stage-A range was frozen before this comparison. Sampled q05/q95 are descriptive quantiles, not probabilities.

## Pre-analysis shape

| Hz | experiment | proxy min | proxy q05 | proxy q50 | proxy q95 | proxy max | status |
|---:|---:|---:|---:|---:|---:|---:|---|
| 10 | 46.1562 | 1.00221 | 1.02487 | 1.26797 | 1.91278 | 2.45364 | outside_sampled_envelope |
| 100 | 4.02433 | 1.00646 | 1.02988 | 1.26533 | 1.88511 | 2.35688 | outside_sampled_envelope |
| 1000 | 1 | 1 | 1 | 1 | 1 | 1 | inside_sampled_q05_q95 |
| 3000 | 0.935658 | 0.50368 | 0.525809 | 0.634503 | 0.832317 | 0.952015 | inside_sampled_min_max |
| 5000 | 0.862552 | 0.405026 | 0.412713 | 0.514176 | 0.709484 | 0.875188 | inside_sampled_min_max |
| 7000 | 0.782724 | 0.346324 | 0.366759 | 0.463719 | 0.694628 | 0.791025 | inside_sampled_min_max |
| 10000 | 0.636483 | 0.300632 | 0.318239 | 0.435804 | 0.678499 | 0.697952 | inside_sampled_q05_q95 |

## Band diagnostics

- **10-100_Hz**: partially_or_fully_outside (2 sampled-anchor points outside min/max).
- **100-1000_Hz**: partially_or_fully_outside (1 sampled-anchor points outside min/max).
- **1-3_kHz**: covered_by_sampled_min_max (0 sampled-anchor points outside min/max).
- **3-10_kHz**: covered_by_sampled_min_max (0 sampled-anchor points outside min/max).

Post-analysis comparison applies the same deterministic 2nd-order 10 kHz Bessel zero-phase ASD factor. It is secondary; detector-side shape interpretation uses pre-analysis.

No parameter range was adjusted and no best-looking member was promoted to a parameter estimate.
