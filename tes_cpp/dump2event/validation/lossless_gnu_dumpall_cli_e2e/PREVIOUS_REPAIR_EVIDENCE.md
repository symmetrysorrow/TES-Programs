# Recovered GNU repair and current run evidence

The earlier serialization repair artifacts were found under the isolated PHITS work area and the repository's private overnight validation data. The repaired-source diff in `patches/dumpall_stream_copy_v2.patch` corrects an initial patch that reused one I/O status variable for both read and write. The final loop uses separate `readstat` and `writestat`, non-advancing 4,096-character reads, `SIZE=`, explicit EOR record completion, and EOF exit.

The following historical artifacts are real and were re-hashed during this audit:

- `phits_WinGfort_lossless_seq.exe`: 44,143,219 bytes; SHA-256 `8bb304d858ba1dd4e94436759e8dc38b7f955a80c301e3c20f6c16bc5709c77e`.
- Its isolated diagnostic source copy `diag_seq_src/ovly15.f`: SHA-256 `1954ef1c1c4f080002eaa1f314428c997c98317ab4fb6d263afc30310f56372e`.
- `diag_seq_src/talls08.f`: SHA-256 `464cc9f5ac274357061b1360121b8a78848337c75e77ec74371322b463188430`.

The independent-seed 1,000-history PHITS run on 2026-10-21 used that executable; the run metadata includes the exact command, executable hash, run input hash, seed, and exit code. It exited 0 and produced a full raw with terminal NCOL=2. The toolchain's detailed source-to-executable build manifest is not preserved, so the source/binary relationship is documented by the isolated directory and hashes, not claimed as independently reproducible.

A new isolated OpenMP build from the installed PHITS source was also attempted, but failed at `mod_ompparallel.o` with `Error: unknown pseudo-op: .tls_common`; this prevents a multi-thread merge claim.
