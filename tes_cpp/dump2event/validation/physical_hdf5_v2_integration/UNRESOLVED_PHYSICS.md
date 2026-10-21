# Unresolved physical scope

The 10,000-history run has 342 unresolved ledger contributions: 333 `unsupported_primary` and 9 `unsupported_secondary`. The 333 primary rows are `ITYP=13` positron records (sampled rows are NCOL=14/JCOLL=14); the nine secondary rows also occur in positron-parent records. Their energies are deliberately invalid/unknown in HDF5. They are not converted to zero.

The HDF5 reconstruction supports only the existing Pb/PHITS 3.341 validated ledger branches. Variable particle weights, full positron stopping and annihilation energy accounting, photon interaction/release branches, all other unrecognized reactions, and parent-child linkage remain outside the validated implementation. The independent 10,000-history histogram has 13 bin differences (L1=18), so no attempt was made to fit or redistribute the missing contributions.

The score representative point for a continuous transport segment is not derived from raw endpoints. It remains absent in HDF5; raw segment start and end are preserved. Parent IDs are also unknown and not guessed.
