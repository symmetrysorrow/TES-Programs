# TES noise mismatch — pulse contamination v4

## Detector v4 statistic

Fixed-noise matched statistic with `C = sigma^2 I`; sigma comes only from pulse pretrigger baselines. Record-local energy normalization is false. Real and injected records share `classify_record_v4()`.

## Pulse-pretrigger noise model

Pulse pretrigger data only; no noise-residual fitting or empirical white/readout floor.

## Raw polarity audit

Measured before orientation. CH0 same-sign fraction is 0.901 and CH1 is 0.972, so the negative-polarity control is excluded from primary validation.

## Block-length provenance

Autocorrelation-derived primary block: 2 ms for CH0 and CH1. Actual 2/5/10 ms sensitivity nulls are stored in `pulse_block_length_selection.json`.

## Primary empirical null

5,000 pulse-pretrigger block-bootstrap records; FWER `1e-3` is resolvable. No periodic tiling.

## Negative-polarity validation

External validation/veto only; not a primary p-value source. Excluded by raw polarity audit.

## Injection recovery

Raw injection followed by exactly one production preprocessing pass. Full pulse q50 recovery is 1.00; full q25 is 0.04; tail q50 at 0–50 ms is 0.26.

## Recovery gate

**False.** `criteria_pass` is computed from measured recovery and null outputs.

## Corrected pulse counts

CH0 accepted 345; definite tail 37, likely tail 26, ambiguous tail 4, pulse-free 278.

## Selection-bias result

Gate: **False**. Correlation is reduced relative to v3 but remains warning-level at 20–100 Hz.

## Pre-analysis all vs pulse-free ASD

All/pulse-free ASD ratios at 10/20/50/100/200 Hz: 0.937/0.927/0.920/0.921/0.941. Pre-analysis has no Bessel.

## Post-analysis ASD

Exactly one 2nd-order 10 kHz Bessel `filtfilt`; simulation transfer uses `|H|^2`.

## PC classification

**PC4** — recovery, selection-bias, and time-reversed-control validation do not satisfy PC3.

## RL classification

**RL5_unvalidated_clean_mask**.

## Auxiliary pulse-free coherence

CH1 is auxiliary/non-production. Pulse-free coherence at 5/10/20/30 Hz is approximately 0.97/0.95/0.82/0.71; phase is near zero.

## 5–30 Hz common-mode evidence

`common_low_frequency_component_supported`: pulse-free data retain the low-frequency common component, without source identification.

## 50–200 Hz local/detector evidence

Coherence falls toward 100 Hz; 50–200 Hz remains a separate local/detector-dynamics diagnostic.

## Common-source topology comparison

**topology_ambiguous**. Bias/readout independent-source transfers are not represented in the current five-state model and exact target parameters remain unresolved. No amplitude fit was performed.

## Clean experiment vs intrinsic TES simulation

Frozen intrinsic-TES comparison is diagnostic only. No noise-residual parameter optimization or physical-source amplitude fit was performed.

## Whether a new physical source may now be added

**No.** Strict target conclusion remains **C — exact target physical case remains unidentified**.
