# Final decision

The explicit Pb EGS5 mode is now reachable from the normal C++ CLI and Python API, while default behavior is unchanged. It validates a narrow configuration, hashes provenance inputs, streams through the existing ledger/HDF5 writer, rejects non-unit weights and unresolved scores, and refuses to overwrite an existing output.

This is **not yet approved for ordinary scientific use**. The key evidence mismatch remains: the earlier parsed raw Dumpall is SHA `d254…`; the later internal diagnostic trace belongs to an October 21 run whose GNU ASCII Dumpall is SHA `c11c…` and is structurally truncated for the current parser. The two later baseline/diagnostic files match each other and their `[t-deposit]` outputs match, but that does not establish a join with the earlier raw file. No new independent seed was run, and no real integrated CLI output was created because the private stopping table is unavailable in this worktree.

Next release gates: regenerate a lossless diagnostic Dumpall from the same exact run as the internal trace; rerun all-addition raw-to-internal matching; provide the local table through an explicit private input; run an independent seed; complete CLI/HDF5 comparison at 100/1,000/10,000 histories; and test the gamma/escape/annihilation conditions. PR #6 remains Draft and unchanged.
