# Energy and position validation

## 100 histories: per-addition internal comparison after HDF5 readback

Readback compared all 38,077 HDF5 deposits against `ALL_TALLY_ADDITIONS.csv`, keyed by the stable addition sequence ID. Matched 38,077/38,077; missing 0; extra 0; event mismatch 0; particle mismatch 0; region mismatch 0; energy mismatch above `1e-12 MeV` 0; maximum energy error `5.7976626971489864e-15 MeV`; positive-position mismatch above `1e-12 cm` 0; maximum checked position error 0 cm. HDF5 sum is `6.929495555108985 MeV`; zero-score rows 130.

The score audit checks raw-derived segment endpoints against the internal trace's interval endpoints and raw secondary points against its validated local-point coordinate. The HDF5 does not label a segment endpoint as its score point. For continuous segments, a separate PHITS representative score coordinate is not populated. Zero additions carry zero energy and do not enter the positive-location score check.

## Independent 1,000 histories

Readback contains 201,971 deposits and 1,000 histories, with 0 unresolved and 0 negative rows. The HDF5-derived per-history energy distribution matches all 150 PHITS `[t-deposit]` bins: L1 difference 0, 58 positive in-range histories on both sides. This is histogram-only verification; internal per-addition and per-history deposit values were unavailable for this run.

## 10,000 histories

Readback contains 1,954,742 tally-ledger rows: 1,954,400 known nonnegative rows and 342 unknown-energy rows. The known-energy HDF5 distribution has 532 positive in-range histories, matching the PHITS histogram's total count, but has 13 differing bins and L1 count difference 18. No exact per-history or per-addition internal comparison is available. Therefore the large run fails distribution reconciliation and cannot support a physics-complete claim.
