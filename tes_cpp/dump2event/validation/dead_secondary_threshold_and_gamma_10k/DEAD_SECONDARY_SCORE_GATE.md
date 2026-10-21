# Dead-secondary score gate

Score rule implemented for the explicitly validated PHITS 3.341 branch: raw child must be dead (`JCLUSTS(4)<0`) and charged (`JCLUSTS(5)!=0`); PHITS reaches the child tally path only if `QCLUSTS(7)>eletb(1)`. Equality is not scored because the source uses strict `.gt.`. The run-specific value is read from `esmin` and `eletb(1)`, not substituted with a generic transport cutoff.

For this run, the input omits `esmin`; PHITS source default and isolated runtime probe both give `esmin=0.001 MeV`, `eletb(1)=9.9999999999999995e-7 MeV`. `emin(12/13/14)=0.001 MeV` are transport cutoffs and are distinct from this tally gate.

The experimental ledger separates:

- `SCORED_RECONSTRUCTED`: PHITS score branch is reached and supported reconstruction applies.
- `NOT_SCORED_BY_PHITS`: source gate proves the scorer is not reached; candidate is emitted only to an audit sidecar, with no deposit row or addition ID.
- `UNRESOLVED`: score status/value/position cannot be established; HDF5 completion is refused.
- `INVALID`: malformed/nonfinite/inconsistent data; processing fails.
- `exact-zero scored addition`: an actual PHITS addition with `tlv=0`; remains a zero-energy deposit row and is not a not-scored candidate.

The threshold API uses the supplied PHITS 3.341 `eletb(1)` value and strict comparison. No clipping, redistribution, or energy normalization is used. The 2 MeV/1000 CLI HDF5 readback is complete: 5,978,912/5,978,912 addition IDs joined, no missing or duplicate IDs, 0 region/position mismatches, 0 unresolved and 1 proven not-scored candidate (NOCAS 398); the 150-bin HDF5 vs PHITS tally L1 difference is 0. The independent-seed 10,000-history run is currently executing; do not treat it as tested until the same-run joins finish.
