# Final decision

HDF5 v2 storage and the experimental streaming path work on three real PHITS 3.341 Dumpalls. The writer transactionality issue is fixed; v2 and v1 Python tests and the existing C++ test are in CI. The 100-history HDF5 readback matches all 38,077 internal additions, energies, regions, and checked positive positions. The independent 1,000-history HDF5 distribution matches the 150-bin PHITS histogram exactly.

The 10,000-history run is partial: 342 records are explicitly unresolved, and the readback histogram differs from PHITS in 13 bins (L1=18). Per-addition internal comparison is unavailable for the 1,000 and 10,000 samples. Parent-child linkage, score points for segments, and full absorption are not established.

**Do not integrate the experimental mode into the default converter or PR #6.** The helper pipeline is currently separate from the normal converter and should remain diagnostic-only. No unresolved or unsupported process is represented as a successful zero-energy deposit. `Pulse_Ms` continues to reject HDF5 v2.
