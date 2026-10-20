# Cutoff and dead-secondary validation

Existing PHITS source audit distinguishes `esmin` from particle-specific `emin(12/13/14)` and the EGS5 internal electron/photon thresholds; the internal dead-secondary score decision is a separate branch. To validate cutoff or dead-child reconstruction, raw fields `JCLUSTS(4)`, `QCLUSTS(7)`, `QCLUSTS(8)`, parent state and scoring region must be parsed and paired to the internal addition.

The GNU files in this run do not pass the first NCOL=4 record, so this turn contributes 0 raw-mapped cutoff/dead-secondary additions and no new reconstruction formula. The existing verified ordinary electron, NCOL=10 boundary, NCOL=11 stop and dead-secondary behavior remains unchanged. Missing internal thresholds must be supplied as explicit input metadata; no guessed value is accepted.
