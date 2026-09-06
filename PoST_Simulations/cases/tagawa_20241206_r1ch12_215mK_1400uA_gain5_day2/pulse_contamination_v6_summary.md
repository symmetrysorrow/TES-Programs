# TES noise mismatch — pulse contamination v6

## Record-wise long-tail conditional null

Each accepted CH0 record uses its own phase-randomized raw-record surrogates; global pooled null is not primary. Sampling is fixed 20→100 with p resolution 0.01.

## Conditional p-value resolution

Option B is selected: primary long-tail alpha is 0.01; 100 surrogates maximum per record.

## Short-axis combined FWER

Contained and edge-aware full-pulse statistics use one combined maximum null; separately calibrated p-value minima are not used.

## Total detector false-positive rate

Short FPR=0.0; long FPR=0.008695652173913044; combined FPR=0.008695652173913044; pass=False.

## Age-resolved injection recovery

Age keys remain separate for 0, 10, 20, 50, 100, and 150 ms; selected-template distributions and lag errors are stored in `pulse_detector_injection_recovery_v6.json`.

## Edge / pre-record recovery

Full-pulse overlap and pre-record tail overlap are stored separately at 25/50/75 percent.

## Recovery gate

Pass=False; measured={"full_q25_full_overlap": 0.0, "full_q50_50pct_overlap": 0.94, "full_q50_75pct_overlap": 0.06, "full_q50_full_overlap": 0.01, "pre_record_tail_age_50ms_q50_50pct": 0.09, "tail_age_0ms_q50": 1.0, "tail_age_100ms_q50": 0.08, "tail_age_10ms_q50": 0.07, "tail_age_20ms_q50": 0.07, "tail_age_50ms_q50": 0.06}.

## Selection-bias result

Short, long, and combined significance-vs-PSD correlations are stored separately; pass=False.

## Corrected pulse/tail counts

CH0 accepted=345; primary exclusive counts={"ambiguous": 11, "full_pulse": 0, "long_tail": 6, "pulse_free_candidate": 328}.

## Whether validated pulse-free exists

No. The detector gate is not passed, so only `pulse_free_candidate` is allowed.

## PC classification

PC4; PC4 is retained while detector validity is unmet.

## Trigger metadata provenance

Target Setting.txt and setting.xml support trigger sample 5000; generic PulseConfig Readout.PreSample=1000 is legacy/conflicting.

## Pulse timing populations

Trigger-relative peak, SNR, amplitude, integral, and half-maximum width are in `pulse_peak_population_audit_v6.json`.

## Pulse template stability

Age-specific source counts, q05/q95 waveforms, and leave-one-out variation are in `pulse_template_stability_v6.json`.

## 5–20 Hz common-mode strength

The two-channel common fraction and principal eigenvectors are diagnostic only; CH1 remains auxiliary.

## Noise principal eigenvector

Saved with 2,000-replicate bootstrap uncertainty.

## Pulse principal eigenvector

Saved using common sample origin and one common scalar normalization per event.

## Bootstrap uncertainty

Magnitude, phase, vector-angle, and common-fraction intervals are saved for both pulse and noise.

## Pulse-vs-noise magnitude agreement

Comparison uses magnitude interval overlap, not cosine similarity alone.

## Pulse-vs-noise phase agreement

Comparison uses circular phase interval overlap and vector-angle uncertainty.

## Final PNS classification

PNS1_pulse_like_common_mode_supported; this does not authorize source-amplitude fitting.

## Common/local spectral decomposition

Largest and second cross-spectral eigenvalues are diagnostic common/local quantities, not physical source fits.

## Independent pulse-event PSD prediction if allowed

Blocked until PNS1 and detector gates pass.

## 50–200 Hz local residual

Partitioned spectra are saved; no intrinsic parameter is changed.

## Clean experiment vs intrinsic TES simulation

Blocked because no validated pulse-free subset exists.

## Whether bath/bias/readout investigation is now justified

No. Adding a stationary physical source is also prohibited. Strict conclusion remains C — exact target physical case remains unidentified.
