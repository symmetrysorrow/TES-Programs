# Remaining limits and provenance cautions

1. Validation covers PHITS 3.341, the stated Pb region/material/medium, unit weights, and the observed EGS5 electron continuous, boundary, cutoff, dead-secondary, and exact-zero additions. It does not validate changing weights, general materials, positron annihilation, photon-local scoring, or all-absorption classification.
2. OpenMP merge was not validated. `USEOMP=true` build failed at `.tls_common`; the `$OMP=2` test used a serial binary.
3. The R1000 HDF5 provenance status intentionally remains `transport_physics_completeness=undetermined`; no TES block allocation or energy-balance claim is made.
4. The normal CLI stdout `physics_config_sha256` differs from the HDF5 `/provenance` resolved-config hash. They appear to hash different configuration layers; this naming mismatch needs a dedicated correction before provenance is treated as canonical. Raw/input/table/HDF5 hashes are separately recorded and agree with the inputs used.
5. The R100 historical run did not preserve its actual executable hash at launch. R1000 has full run metadata, including executable hash.
6. The 1000-history PHITS output histogram match is same-run and not an independent statistical comparison. A 10,000-history normal CLI HDF5 run was not performed.
7. The recovered HATCH table matches the known exact SHA and medium-2 Pb rows, but the original HATCH generation transcript is historical and was not rerun.
