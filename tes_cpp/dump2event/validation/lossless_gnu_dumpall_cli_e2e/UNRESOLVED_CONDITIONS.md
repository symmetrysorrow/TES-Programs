# Unresolved conditions

1. Recover the exact PHITS writer patch/source and matching binary from the prior 100-history “lossless” experiment, or rebuild and validate a new isolated implementation.
2. Preserve record contents and requested EOF/newline semantics; the standalone probe did not preserve an unterminated final line exactly.
3. Verify actual thread-file merge ordering, complete NCOL=4 tuples, all record boundaries, and NCOL=2 termination.
4. Recover and hash-verify the private Pb HATCH table (`20fd86eb…514e9`) or regenerate it from explicit PHITS 3.341 HATCH settings and validate coefficient values.
5. Recover/recreate the diagnostic source hooks and executable so internal additions and raw Dumpall originate from the same run.
6. Run the staged 100 / 1,000 / 10,000-history sequence, normal CLI HDF5 v2 output, HDF5 readback, addition/position/region matching, histogram comparison, and independent seed.

No PHITS source, installation executable, original Dumpall/input, or existing HDF5 was modified. No unsupported field was imputed.
