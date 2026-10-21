# Negative, zero, and duplicate addition audit

- R100: 38,077 unique HDF5 addition IDs; 130 exact-zero rows; 0 negative; 0 unresolved.
- R1000: 195,185 unique HDF5 addition IDs; 1,141 exact-zero rows; 0 negative; 0 unresolved.
- PHITS `TALLYSEQ` contains 30 repeated record serials in R100 and 126 in R1000 because a single raw record can lead to multiple scorer increments. They were matched by occurrence order; no HDF5 addition was missing, extra, or overmatched.
- `output=deposit` excludes zero-total histories from the nonzero histogram distribution: 91/100 and 944/1,000 histories had total deposit exactly zero; positive histories were 9 and 56.

Zero energy rows remain explicit in HDF5. No clip, forced normalization, or synthetic nonzero deposit was used.
