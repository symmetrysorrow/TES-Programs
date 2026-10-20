# Weighted addition validation

PHITS 3.341 internal diagnostics were obtained for fixed source weights 0.5, 1 and 2, each at 100 histories. Each run reports 22,069 internal score calls and an unweighted sum of `tlv = 3.7357872960843292 MeV` (21,955 positive, 114 zero, zero negative). Existing normalized 150-bin histograms coincide across those three cases. The PHITS source shows the unweighted `tlv` accumulated within the history and history-end histogram scoring uses `sweight * ratio`; a warning is emitted if OLDWT changes within a history.

These are internal diagnostics and history-level weighted behavior, not a raw reconstruction validation. The GNU Dumpall files fail at their first NCOL=4 particle tuple, so the requested raw OLDWT/WT/QCLUSTS(8) to `sweight`/`ratio`/`tlv` per-addition mapping has **0 verified records** in this run. No fixed-weight or changing-weight branch was newly implemented. Prior results remain evidence about the PHITS side only.
