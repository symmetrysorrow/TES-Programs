# Lossless Dumpall and terminal-record checks

The serial diagnostic runs used the GNU Fortran non-advancing record-copy repair and the exact-sequence diagnostic build. The repair preserves logical Dumpall record payloads without numeric inference; it writes to run-local output only. The 100-history samples have successful PHITS exit code 0, 100 NCOL=4 history headers, and terminal NCOL=2; the 0.511 MeV sample has 250,905 parsed records and the 2 MeV sample has 549,674. The raw-record sequence sidecar maps every track record used by the same-run join.

The 1,000-history 0.511 MeV run exited 0 and produced 2,541,373 parsed records across 1,000 histories with terminal NCOL=2. Its Dumpall is 3,594,947,832 bytes. The normal HDF5 CLI reported the same record/history counts and `NCOL=2` terminator. Its tally and addition-level comparison are still in progress. The 2 MeV 1,000-history run has not yet returned at report update time.

The 10,000-history Co-60 baseline is not accepted by the converter: the strict parser reports line 5, NCOL=4, particle-header field with 11 numeric values where 8 are required. The original was not rewritten. This is an input serialization/parser failure, not evidence about tally physics. No truncation, field guessing, or HDF5 publication was performed.
