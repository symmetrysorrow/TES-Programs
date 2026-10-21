# GNU `cpcasdump` record truncation

Audit date: 2026-10-21 (Asia/Tokyo). The installed PHITS source `C:\PHITS\src\ovly15.f` was read only; its original file SHA-256 was `573038237a089c40ac4a3ab53ab5e9ce24345ef58959fe9c06132f6d85152bb8`.

The GNU ASCII pass-through in `cpcasdump` used `character(80) dummy`, then formatted `(a)` input/output. The field loss is in PHITS serialization, before the converter sees the data. The repaired source copy uses a 4,096-character **chunk** and non-advancing reads with `SIZE=`. It loops while `iostat=0`, writes each transferred chunk unchanged, emits a record terminator only on `iostat_eor`, and exits at `iostat_end`. Separate read and write status variables prevent the earlier status-clobbering defect. No numeric token is parsed or reformatted.

Evidence:

- Repaired 100-history raw: 647,830 physical lines; longest line 338 characters; 267,683 lines exceed 80; none is exactly 80; final line is `2` (NCOL=2). C++ parser: 38,345 structural records, 100 histories, terminal NCOL=2.
- New independent-seed 1,000-history raw: 3,332,756 lines; longest line 338; 1,378,712 lines exceed 80; none exactly 80; final line is `2`. C++ parser: 198,060 records, 1,000 histories, terminal NCOL=2.
- Synthetic GNU Fortran 16.1.0 probe copied CRLF input containing 6 logical records of lengths 5, 80, 81, 12,000, 0, and 22 characters. All record payloads and boundaries matched. The unterminated final record gained the runtime's final record terminator, so byte-for-byte newline preservation is **not** claimed.
- `C:\PHITS` source was not changed.

Source artifact from the earlier isolated PHITS build: `diag_seq_src/ovly15.f` SHA-256 `1954ef1c1c4f080002eaa1f314428c997c98317ab4fb6d263afc30310f56372e`; the used diagnostic executable SHA-256 is recorded in `RUN_IDENTITY_AND_SHA_AUDIT.md`. A fresh OpenMP-enabled build from an isolated copy of `C:\PHITS\src` was attempted on 2026-10-21 but failed in `mod_ompparallel.o`: GNU assembler rejected `.tls_common`. Therefore multi-thread cpcasdump merge order is **not validated**. A `$OMP=2` run on the non-OpenMP executable is not counted as a threaded test.

See `LOSSLESS_WRITER_IMPLEMENTATION.md`, `THREAD_MERGE_VALIDATION.md`, and `PREVIOUS_REPAIR_EVIDENCE.md`.
