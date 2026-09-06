# TES noise mismatch — pulse contamination v5

## v5 preprocessing correctness
Raw bootstrap and processed covariance domains are separated; each processed path receives one production preprocessing pass.

## Trigger-relative pulse timing
Target `Setting.txt` establishes trigger sample 5000 and peak window `[5000, 15000)`. Legacy `PulseConfig.PreSample=1000` is recorded as conflicting metadata.

## Raw polarity result
See `raw_pulse_polarity_audit_v5.json`; polarity is evaluated only in the trigger-compatible window.

## Tail-age template bank
Aged templates 0/10/20/50/100/150 ms are stored in `pulse_template_library_v5.json`.

## Edge detector validation
Available overlap and template-only denominator are used with minimum overlap 0.25.

## Short-timescale null
Short/full pulses use the raw pretrigger bootstrap null.

## Long-tail conditional surrogate null
Long tails use 999 phase-randomized, simulation-blind surrogates preserving per-record PSD magnitude.

## IAAFT validation
Secondary IAAFT validation is stored separately and does not set thresholds.

## Injection recovery / Recovery gate
Recovery gate: **False**.

## False-positive controls / Selection bias
Time-reversed controls remain frozen-threshold validation; selection-bias pass: **False**.

## Corrected pulse/tail counts
CH0 accepted `345`; counts are in `noise_record_pulse_classification_v5.json`.

## All vs pulse-free candidate ASD
No `validated_pulse_free` subset exists because the detector gate is not passed.

## Final PC classification
**PC4**.

## 5–30 Hz auxiliary common mode / cross-spectral principal eigenmode
CH1 remains auxiliary/non-production. Pulse-free candidate common coherence is retained independently of detector validity.

## Paired pulse channel signature / Pulse vs noise common-mode similarity
Comparison status: **pulse_channel_signature_match**; no amplitude fit.

## 50–200 Hz intermediate-pole diagnostic
Existing pulse-only timescales are reported without fitting noise residual poles.

## Whether common forcing is pulse-like
Inconclusive until detector validity and signature separation pass.

## Whether bath/bias/readout topology investigation is justified
Not yet; exact target topology remains unresolved.

## Clean experiment vs intrinsic TES simulation
Blocked because validated pulse-free does not exist.

## Whether adding a new stationary source is allowed
**No.** Strict target conclusion remains **C — exact target physical case remains unidentified**.
