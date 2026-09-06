# FINAL PULSE CLOSURE AUDIT

The full fixed-M recomputation was not completed within the available compute window. No incomplete result is represented as PASS. Existing v6 evidence already has detector validity FAIL, selection-bias FAIL, PC4, and no validated pulse-free subset.

## Final long-tail conditional detector

NOT RUN to completion. The final implementation specifies fixed M=199, record-local phase-randomized nulls, and no optional stopping.

## Final short-mode calibration

NOT RUN to completion. The implementation specifies separate contained/left-edge/right-edge calibration and an independent validation combined null.

## Final control FPR

NOT RUN to completion; sign-inverted control was not completed. The final policy is short alpha 0.001, long alpha 0.01, union upper scale 0.011, and 95% binomial upper bound below 0.02.

## Final injection recovery

NOT RUN to completion. The implementation specifies an own fixed-M conditional null per injected record and no pooled injection null. Existing v6 recovery already fails grossly (full q50 100% = 0.01 versus full q50 50% = 0.94).

## Final detector validity

FAIL / unavailable for final fixed-M run. `validated_pulse_free` is not permitted.

## Final selection bias

FAIL in existing v6 evidence: final combined score passed its combined diagnostic, but the required final audit was not completed and the short-only correlation remained high.

## Validated pulse-free availability

No.

## All vs clean removed power fraction

UNAVAILABLE; no validated clean subset and no completed final PSD bootstrap.

## Pulse/noise topology equivalence

PNS3. CI overlap alone is not equivalence, and independent margins/direct-difference equivalence were not completed.

## Coherent CH0 power

The final artifact defines `|S01|²/S11` and `S00-|S01|²/S11`; no source amplitude fit is made.

## Pulse-event predicted PSD

UNAVAILABLE because detector validity failed.

## Pulse predicted / measured common power

UNAVAILABLE.

## Final PC classification

PC4.

## Final PNS classification

PNS3.

## FINAL PULSE DECISION

FINAL PULSE DECISION:

CLOSE_PULSE_TRACK_UNRESOLVED_WITH_EXISTING_DATA

Should another pulse detector version be built?

NO

Pulse may exist, but existing data cannot support a quantitative pulse-dominance determination. Further same-data detector iteration is not scientifically justified. Proceed to amplitude-free 5–30 Hz common-mode topology investigation. Strict target conclusion remains C — exact target physical case remains unidentified.
