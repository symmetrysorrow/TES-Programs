# Pulse contamination summary

## Pulse contamination rate

CH0: `101` definite, `5` ambiguous, `239` pulse-free of 345 accepted records; definite fraction `29.275%`.
CH1 has only 1 independently accepted record under the unchanged production acceptance predicate.

## Detection confidence / false positives

Primary FPR is fixed at 1e-3 from pulse pretrigger null portions. Forward-control false-positive rates: CH0 `0.001667`, CH1 `0.001667`. Time-reversed and sign-inverted controls are retained in `pulse_false_positive_audit.json`.

## Pulse-free record count

CH0 strict clean count: `239`. No pulse subtraction is used.

## All vs pulse-free ASD

CH0 ASD all/clean ratios at anchors: `{"10": 1.0106473230645288, "100": 1.0090002520248031, "1000": 0.9643156039350326, "10000": 1.016446635408553, "20": 0.9814401046196577, "200": 1.017465436934332, "3000": 0.996817209303773, "50": 0.9663565680996162, "500": 0.9892738215769359, "5000": 1.0195086899879964, "7000": 0.9874088976631673}`. The low-frequency change is small.

## Contaminated subset ASD

Reconstructed pulse-only power fractions of all-record power at 10/20/50/100/200 Hz: `{"10": 0.00029838820267820225, "100": 0.013951400551316237, "20": 0.001369314341721452, "200": 0.044485333901206306, "50": 0.0064003056196037175}`.

## Coherence before/after pulse rejection

Exact-key paired coherence is inconclusive because independent acceptance leaves `1` pair(s).

## CH0/CH1 pulse coincidence

Detected paired noise pulses: `0`; pulse-dataset amplitude-ratio and lag diagnostics are in `pulse_channel_coincidence.json`.

## Predicted pulse-only PSD

The pulse-only reconstruction uses detected event lag, template class, and time-domain amplitude only. No spectral residual scaling is fitted.

## Stationary simulation vs pulse-free experiment

Pulse-free CH0 vs pulse-gated model covers `6/11` anchors; the 10–100 Hz excess remains in the clean subset.

## Stationary + pulse contamination vs original experiment

Additive predicted/all relative errors at anchors: `{"10": -0.010384380766577772, "100": -0.001906303835016221, "1000": 0.043897363946795, "10000": -0.016145831256855203, "20": 0.019582607905155758, "200": 0.005210788570430669, "3000": 0.0038880092820328116, "50": 0.03790260008740143, "500": 0.02112074309025802, "5000": -0.01891978123334881, "7000": 0.01282698245621483}`. Agreement is a consequence of the small pulse contribution.

## Remaining residual

The remaining 10–100 Hz excess is not explained by detected pulse contamination. CH0/CH1 coherence classification remains inconclusive because of the strict independent CH1 acceptance count.

## Final PC classification

**PC3 — pulse contamination negligible** for the observed CH0 low-frequency excess, conditional on the fixed-FPR detector. The previous R5 provisional label is superseded.

Strict target conclusion: **C — exact target physical case remains unidentified**.
