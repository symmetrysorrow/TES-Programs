# TES noise mismatch — pulse contamination v3

## Detector v3 correctness

Raw and one-pass processed pulse templates are stored separately. The detector uses full-record, all-lag, all-template statistics and conservative dual-null family p-values. The `1e-4` category is removed.

## Raw-domain injection recovery

Raw pulse is injected into raw bootstrap baseline, then production preprocessing is applied exactly once. Recovery criteria pass: **False**.

## Empirical FWER resolution

Target primary FWER is `1e-3`; operational rank resolution is recorded in `pulse_detection_decision_policy_v3.json`. Negative-polarity control has its actual finite-record resolution.

## Block-bootstrap null

Pulse-only 5 ms block-bootstrap null; no periodic baseline repetition construction.

## Negative-polarity null

Accepted noise records scanned with physically impossible negative-polarity templates; it is not a simulation or spectrum fit.

## Dual-null consistency

Classification uses the conservative maximum of block-bootstrap and negative-control family p-values.

## Full/tail/edge recovery

See `pulse_detector_injection_recovery_v3.json`; the recovery gate is not passed, so the clean mask is not definitive.

## Corrected pulse counts

CH0 accepted `345`; counts: `{"ambiguous_full": 1, "ambiguous_tail": 50, "definite_full": 0, "definite_tail": 0, "likely_full": 1, "likely_tail": 0, "pulse_free": 293}`.

## Selection bias

Spearman correlations between `−log10(p)` and `log10(PSD)` are stored in `pulse_selection_bias_audit.json`; the measured 10–100 Hz correlations are materially negative (approximately −0.46 to −0.54), so the clean mask is selection-sensitive. No selection feedback is used.

## Pre-analysis all vs clean ASD

Stored explicitly in `pulse_partitioned_noise_spectra_v3.json`; pre-analysis has no Bessel filter. Clean/all ASD ratios at 10/20/50/100/200 Hz are `{"10": 0.9394, "100": 0.9535, "20": 0.9348, "200": 0.9839, "50": 0.9382}`.

## Post-analysis all vs clean ASD

Stored separately with one 10 kHz Bessel `filtfilt` pass.

## 10–200 Hz slope

The clean-v3 candidate pre-analysis slope over 10–200 Hz is `γ=-1.0546`. It remains diagnostic only until recovery and mask-bias criteria pass.

## Record-length scaling

Not promoted: `record_length_scaling_v3.json` is blocked because the clean mask is not recovery-valid.

## Detrend dependence

Mean-removal, linear, and quadratic diagnostics are stored; they do not replace production spectra.

## Auxiliary CH0/CH1 cross-spectrum

Stored as auxiliary only. CH1 is not production-accepted; exact-key hardware-valid pairs are kept separate.

## Pulse-free coherence

Auxiliary pulse-free coherence is reported without calling it production-accepted coherence.

## Clean vs TES simulation

No PC3 promotion or stationary-source addition is allowed before recovery and mask-bias criteria pass.

## Final PC classification

**PC4** — detector recovery/null/mask validity is insufficient.

## Final RL classification

**RL5_unvalidated_clean_mask** — record-length classification is not promoted.

## Whether stationary physical-source investigation is now allowed

**No.** PC3 + RL1 conditions are not met. Strict target conclusion remains **C — exact target physical case remains unidentified**.
