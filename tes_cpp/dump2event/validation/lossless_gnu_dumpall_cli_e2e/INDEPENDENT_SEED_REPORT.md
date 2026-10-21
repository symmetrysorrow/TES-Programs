# Independent-seed R1000 validation

A new PHITS 3.341 run used a different 64-bit seed from the earlier R100 run while keeping the Pb geometry/material, source lines, transport settings, cutoff, tally, and medium-2 stopping table fixed. `maxcas=1000`; PHITS exit code 0. The run's exact input, seed, executable hash, raw hash, and HDF5 hash are in `RUN_IDENTITY_AND_SHA_AUDIT.md`.

The ordinary `dump2event --physics pb-egs5 --hdf5-version 2` conversion matched all 195,185 same-run PHITS internal additions. Maximum absolute energy error was `5.830839283627043e-15 MeV`; per-event max sum error was `5.129230373768223e-14 MeV`; region and position mismatch counts were 0. The `[t-deposit]` 150-bin distribution (0.01 MeV bins, `[low, high)`) matched exactly after excluding zero-deposit histories, which the PHITS positive-score histogram does not include: 56 positive histories on each side, L1 count difference 0, max per-bin difference 0, underflow 0, overflow 0.

This is a same-run cross-reconstruction check, not an independent statistical test of two PHITS runs. The PHITS tally and raw/HDF5 artifacts share the same histories. It validates the observed seed and the scoped physics branches only.
