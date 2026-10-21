# Lossless writer implementation and observed limits

The accepted serialization repair is the non-advancing chunk loop in `cpcasdump` described in `GNU_CPCASDUMP_ROOT_CAUSE.md`. The 4,096-character buffer is a repeated I/O chunk, not a maximum record length; records longer than the chunk are read in multiple iterations. It preserves logical record content without numeric parsing.

## Executed validations

- GNU Fortran 16.1.0 synthetic copier: 5-, 80-, 81-, 12,000-character records, empty record, and unterminated final record. All 6 record payloads and boundaries matched; the final record received a newline on output.
- Repaired PHITS 3.341 100-history raw: 38,345 records / 100 NCOL=4 headers; terminal NCOL=2; max line 338; 267,683 lines above 80.
- Repaired PHITS 3.341 independent-seed 1,000-history raw: 198,060 records / 1,000 headers; terminal NCOL=2; max line 338; 1,378,712 lines above 80.
- 100-history no-hook and diagnostic-sequence runs with the same seed and writer have identical Dumpall SHA-256 and `[t-deposit]` SHA-256. This supports non-perturbation for that tested serial configuration.

## Not established

- Exact byte preservation of CRLF versus LF or of an absent final newline; Fortran formatted record I/O normalizes the final terminator.
- OpenMP/threaded merge order. The installed makefile defaults to `USEOMP` off. A fresh `USEOMP=true` build failed at assembler `.tls_common`; no threaded PHITS output is claimed.
- A fresh rebuild of the exact historical diagnostic executable from the recorded source tree; the new 1,000-history run records its binary hash, but not a new source-to-binary reproducible build log.
- 10,000-history conversion through the normal HDF5 v2 CLI in this turn.
