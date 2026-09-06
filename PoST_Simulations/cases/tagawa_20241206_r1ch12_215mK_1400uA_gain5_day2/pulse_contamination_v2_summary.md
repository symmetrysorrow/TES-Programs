# TES noise mismatch — pulse contamination v2

## Detector v2 validity

Record-wise empirical null calibration, full-rate tail templates, all-lag scanning, and family-wise template significance are implemented. Injection recovery is insufficient for an unconditional clean mask; the result is therefore PC4, not frozen PC3.

## Record-wise false positive calibration

Nulls use pulse pretrigger baselines extended to 0.2 s and the production preprocessing. The primary family-wise target is FWER `1e-3`; calibration used `1024` CH0 null records and includes all templates and valid lags.

## Tail-only recovery efficiency

CH0 tail-only efficiencies by tail age (ms): `{"0": 0.08333333333333333, "10": 0.08333333333333333, "100": 0.08333333333333333, "150": 0.08333333333333333, "50": 0.16666666666666666}`. Full-pulse efficiency is `0.167`, and onset-near-end efficiency is `0`.

## Corrected pulse contamination fraction

CH0 accepted records: `345`; v2 counts: `{"ambiguous": 6, "definite_full_pulse": 0, "definite_tail_only": 0, "likely_pulse": 176, "pulse_free": 163}`. The 176 likely records are not promoted to a definitive contamination fraction because recovery is incomplete.

## All vs clean_v2 ASD

Clean_v2 contains `163` records. Clean ASD anchors: `{"10": 0.00096214949588305, "100": 8.28434457682301e-05, "1000": 1.655820792288282e-05, "10000": 3.2257594630528893e-06, "150": 5.37745719633772e-05, "20": 0.0004332492781768366, "200": 3.835132653399284e-05, "30": 0.00029442016954407105, "300": 2.391665185232059e-05, "3000": 1.3747350214435986e-05, "5": 0.0012496135925359362, "50": 0.00017277309494244888, "500": 2.4398187692080298e-05, "5000": 1.0810569940042045e-05, "70": 0.00011687873726302177, "7000": 7.520838296982709e-06}`. No pulse subtraction is used.

## 10–200 Hz spectral slope

Descriptive clean_v2 gamma: `-1.0644`; this is close to ASD proportional to 1/f, but is not a source fit.

## Record-length scaling

The maximum common-bin ratio between 0.2/0.1/0.05 s is `1.252`; conditional classification: **RL1_conditional_on_clean_v2** (1.5 stability cutoff).

## Detrend dependence

Linear detrending is negligible over 10–200 Hz; quadratic detrending changes the lowest-frequency anchor more than the mid-band. Detrend outputs are diagnostic only and do not replace production ASD.

## Per-record low-frequency dominance

Top-1% records contribute approximately 7–9% of 10/20/50 Hz power; this is not a single-record-dominated excess.

## CH1 acceptance provenance

CH1 remains 1 production-accepted record. No CH1-specific dynamic-range, ADC-scale, gain, or saturation rule was found; no production-accepted coherence is claimed.

## Corrected coherence if available

Unavailable as production-accepted coherence because CH1 acceptance provenance is unresolved; the auxiliary paired-channel audit is kept separate.

## Clean_v2 vs TES simulation

The existing frozen pulse-consistent simulation ensemble was compared descriptively; no parameter generation, residual fit, or scenario promotion was performed. The 10–200 Hz clean excess remains outside the frozen envelope.

## Final PC classification

**PC4** — detector uncertainty / long-tail ambiguity is material; PC3 cannot be retained as a definitive result.

## Final RL classification

**RL1_conditional_on_clean_v2** — the record-length result is conditional on the v2 clean mask; it is not evidence that the physical source is identified.

## Next physical hypothesis

Do not advance to bias/bath/readout source attribution yet. First improve tail-only and end-edge injection recovery using pulse-only morphology and a defensible record-wise null. Strict target conclusion remains **C — exact target physical case remains unidentified**.
