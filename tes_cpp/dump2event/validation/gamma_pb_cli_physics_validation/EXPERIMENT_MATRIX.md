# Pb mono-gamma and HDF5 v2 validation matrix

Run scope: PHITS 3.341, Pb region 1000/material 2, EGS5 medium 2, unit source and particle weights, electron/positron/photon cutoff 1 keV, `[t-deposit] mesh=reg, reg=1000, part=all, deposit=0, unit=3, output=deposit, e-type=2, 150 bins from 0 to 1.5 MeV`. The mono-gamma inputs are derived copies of the earlier isolated templates; original PHITS input and installation files were not edited.

| Beam / histories | PHITS run | ordinary CLI HDF5 v2 | Internal addition join | Tally result |
|---|---|---|---|---|
| 0.511 MeV / 100 | completed, exit 0 | completed, exit 0 | 250,938 / 250,938; zero identity, region and position mismatches | 150-bin counts exact; 2 zero-deposit histories excluded; no overflow |
| 2 MeV / 100 | completed, exit 0 | completed, exit 0 | 549,639 / 549,639; zero identity, region and position mismatches | 150-bin counts exact for in-range events; 9 zero-deposit and 19 >1.5 MeV overflow histories |
| 0.511 MeV / 1,000 | completed, exit 0; 2,541,373 records and terminal NCOL=2 | completed: 1,000 histories, 2,541,564 deposits, no negative/unresolved | all 2,541,564 additions joined; zero event/particle/region/position mismatches | exact 150-bin match; 27 zero histories |
| 2 MeV / 1,000 | completed, exit 0; PHITS reports 1,000 cases and zero remaining batches | CLI failed closed (exit 1): 5,978,912 reconstructed rows plus 1 below-guard unresolved candidate; no HDF5 published | raw-ledger ordered comparison: 5,978,912 / 5,978,912 matched; zero identity mismatches; max energy diff 4.68e-14 MeV | raw-ledger event histogram exact: 150-bin L1=0; 87 zero and 275 overflow histories |
| Co-60 / 10,000 | existing completed raw was tried as a CLI input; failed closed at malformed NCOL=4 layout (see `TEN_THOUSAND_HISTORY_CLI_REPORT.md`) | failed, no HDF5 published | not applicable | existing tally not used to claim a new gamma-run comparison |

The 100-history comparison uses the same PHITS run's PRE/POST tally trace, `TALLYSEQ`, lossless raw-record sequence sidecar, repaired ASCII Dumpall, and HDF5 v2. The histogram comparison is a same-run distribution comparison; it does not establish individual-event physical completeness beyond the tested slice. Details and hashes are in the linked reports.
