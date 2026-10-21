# Final decision

**The normal CLI success objective is demonstrated for the validated 100- and 1,000-history PHITS 3.341 unit-weight Pb EGS5-electron slice.** Both produced readable HDF5 v2; all same-run additions matched the PHITS internal trace within about `6e-15 MeV`, region and recorded position mismatches were zero, and the positive-history `[t-deposit]` histogram matched all 150 bins exactly.

This is not approval for general or default production use. The current HDF5 provenance keeps physics completeness and absorption status undetermined. Variable weights, positron/photon scoring, other materials, OpenMP merge, and 10,000-history CLI conversion remain unverified. The resolved-config hash naming mismatch also needs cleanup. Keep `pb-egs5` explicit and experimental; do not fold this into PR #6 or the default converter based solely on this scoped result.
