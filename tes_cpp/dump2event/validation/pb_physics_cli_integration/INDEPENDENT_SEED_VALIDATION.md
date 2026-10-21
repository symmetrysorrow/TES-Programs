# Independent seed validation

**Not completed.** No new independent seed was chosen and then run in both a normal and non-perturbative diagnostic PHITS 3.341 build during this integration task. The October 21 baseline/diagnostic pair reused the earlier 64-bit seed; it is a diagnostic non-perturbation pair, not an independent-seed validation.

The historical same-seed pair has identical `[t-deposit]` SHA-256 and identical Dumpall SHA-256 within that pair. However, the GNU diagnostic Dumpall is truncated at the first NCOL=4 particle tuple under the current structured reader, so it cannot support same-run raw-only score reconstruction. A repaired Dumpall and the local HATCH table are prerequisites for the requested same-run comparison.

No positron-branch count, per-addition error, region/position agreement, or histogram L1 is claimed for an independent seed. New-seed physical reproducibility remains a release gate.
