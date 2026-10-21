# OpenMP status

OpenMP merge validation was outside this run matrix. The previously recorded isolated PHITS GNU OpenMP build failed at assembler pseudo-op `.tls_common`; the current 100- and 1,000-history PHITS runs use the serial exact-sequence executable with `OMP_NUM_THREADS=1`. No threaded PHITS behavior or merge-order equivalence is claimed.
