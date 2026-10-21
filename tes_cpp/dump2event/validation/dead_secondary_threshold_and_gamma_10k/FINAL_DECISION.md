# Interim implementation decision

Scope currently established: PHITS 3.341, 2 MeV mono-gamma, Pb region 1000, unit weight, 1 keV particle cutoffs, and the explicit experimental `pb-egs5` HDF5 v2 mode.

PHITS source plus same-run runtime measurement establish the dead-charged-secondary gate `QCLUSTS(7) > eletb(1)`, with default `esmin=0.001 MeV` and `eletb(1)=9.9999999999999995e-7 MeV` for this run. The below-threshold NOCAS 398 child is proven not scored and is recorded in an audit sidecar rather than fabricated as a zero deposit. The standard PHITS internal addition at that parent reaction is separately retained.

The 2 MeV / 1000-history ordinary CLI HDF5 v2 output matched all 5,978,912 same-run internal additions; no ID, event, particle, region, or coordinate mismatches; maximum energy error `4.680023729664029e-14 MeV`; 150-bin tally L1 difference 0. This supports this bounded experimental slice only.

The requested independent-seed 10,000-history gamma run was stopped at 2026-10-22 01:34 JST by terminating only verified PHITS PID 41452, as requested. The driver recorded exit code -1; the partial files are not a normal PHITS run. Its normal completion, Dumpall terminator, HDF5 v2 output, same-run per-addition correspondence, and history distribution remain unverified. The bounded 2 MeV/1,000-history dead-secondary case was re-audited separately; see `STOPPED_10K_AND_1000_REREAD_20261022.md`. This is **not a final pass**, not a production-ready physical engine, and not grounds to change PR #6's default deposit formula.
