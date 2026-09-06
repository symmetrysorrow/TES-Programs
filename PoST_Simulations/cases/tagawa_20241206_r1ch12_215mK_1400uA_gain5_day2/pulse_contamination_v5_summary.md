# TES noise mismatch — pulse contamination v5

## v5 preprocessing correctness
Raw bootstrap and processed covariance domains are separate; the parity artifact reports one production preprocessing pass and no double filtering.

## Trigger-relative pulse timing
Trigger is sample 5000 and the independently justified peak window is samples 5000–15000. Whole-record absolute-maximum peak search is not used.

## Raw polarity result
CH0/CH1 polarity is evaluated on trigger-window signed peaks and integrals; see `raw_pulse_polarity_audit_v5.json`.

## Tail-age template bank
Measured age-specific templates exist for 0, 10, 20, 50, 100, and 150 ms in both raw and processed domains.

## Edge detector validation
The detector supports fully contained, onset-near-end, and pre-record negative-lag overlap with template-only available-overlap normalization and a fixed minimum overlap of 0.25.

## Short-timescale null
Full/short detection uses the 5,000-record raw pretrigger bootstrap null.

## Long-tail conditional surrogate null
Long-tail detection uses 1,000 simulation-blind Fourier phase-randomized surrogates; per-record non-DC PSD magnitude is preserved within the recorded tolerance.

## IAAFT validation
100 IAAFT surrogates per channel are secondary validation only and do not set thresholds.

## Injection recovery
Injection is raw-domain, followed by exactly one production preprocessing pass, and all decisions use `classify_record_v5()`.

## Recovery gate
`criteria_pass = False`; measured values: `{"false_positive_per_record": 0.001, "full_q25": 0.03, "full_q50": 1.0, "onset_near_end_q50": 0.01, "pre_record_tail_q50": 0.0, "primary_threshold_rho": 224.9168805572623, "tail_q25_age_0_to_50ms": 0.0, "tail_q50_age_0_to_50ms": 0.16666666666666666, "tail_q50_age_100ms": 0.16666666666666666}`.

## False-positive controls
Time-reversed controls use frozen primary thresholds; status is `pass`.

## Selection bias
Pass is `False`; 10–100 Hz correlations remain above the predeclared |rho|<0.2 gate.

## Corrected pulse/tail counts
CH0 accepted `345`: `{"long_ambiguous": 15, "long_definite": 0, "long_likely": 1, "primary_ambiguous": 16, "primary_full_pulse": 0, "primary_long_tail": 0, "pulse_free_candidate": 329, "short_ambiguous": 2, "short_definite": 0, "short_likely": 0}`. Primary subsets are exclusive.

## All vs pulse-free candidate ASD
Candidate count is `329`; all/candidate ASD ratios are `{"post_analysis": {"10": 1.0737422869548812, "100": 1.0721373009980906, "1000": 1.0068572035774583, "10000": 1.0206525569873404, "150": 1.0962037216679252, "20": 1.062362452255824, "200": 1.0951853323686258, "30": 1.0667983151684388, "300": 1.069077500907365, "3000": 0.9896496806822905, "5": 1.0717556046943528, "50": 1.0532036877699948, "500": 1.0026114145246532, "5000": 0.9985960143731735, "70": 1.0642430814158756, "7000": 0.9959926248926886}, "pre_analysis": {"10": 1.0737422787882607, "100": 1.0721372940450118, "1000": 1.0068571378535398, "10000": 1.0206644162766771, "150": 1.0962036507253412, "20": 1.0623624700298626, "200": 1.0951849764594153, "30": 1.0667983108078969, "300": 1.0690770695860035, "3000": 0.989647520279577, "5": 1.0717556128537864, "50": 1.0532037318082474, "500": 1.0026106829796944, "5000": 0.998599284258496, "70": 1.0642431985028964, "7000": 0.9959930533171689}}`.

## Whether validated pulse-free exists
No: `validated_pulse_free` is not created before the detector gate passes.

## Final PC classification
`PC4`; PC4 remains because the detector recovery gate is not passed.

## 5–30 Hz auxiliary common mode
Exact-key CH0/CH1 auxiliary coherence is `{"10": 0.9459061709030585, "100": 0.09537244529099541, "20": 0.8034145629963291, "30": 0.6835186051773136, "5": 0.9681962942265714, "50": 0.39218755566164154}` and remains strong in the candidate subset.

## Cross-spectral principal eigenmode
Saved in `noise_common_mode_eigenanalysis.json`; CH1 is explicitly auxiliary/non-production.

## Paired pulse channel signature
Saved from exact-key trigger-relative paired pulse windows in `paired_pulse_channel_signature.json`, without noise PSD input.

## Pulse vs noise common-mode similarity
Classification is `pulse_channel_signature_match` with global-phase-invariant similarity values `{"10": 0.9382121189166528, "100": 0.9051969885300822, "20": 0.9396385050912888, "30": 0.9516570880749495, "5": 0.98438674777311, "50": 0.8096353430664094, "70": 0.9330910527459075}`; no amplitude fit.

## 50–200 Hz intermediate-pole diagnostic
Existing pulse-only timescales are compared with the band; no residual pole fit is performed.

## Whether common forcing is pulse-like
Not established: signature comparison is descriptive, while detector validity and selection-bias gates fail.

## Whether bath/bias/readout topology investigation is justified
Not yet. Exact target physical case and independent topology linkage remain unresolved.

## Clean experiment vs intrinsic TES simulation
Simulation comparison is blocked because a validated pulse-free subset does not exist.

## Whether adding a new stationary source is allowed
No. Strict target conclusion remains **C — exact target physical case remains unidentified**.
