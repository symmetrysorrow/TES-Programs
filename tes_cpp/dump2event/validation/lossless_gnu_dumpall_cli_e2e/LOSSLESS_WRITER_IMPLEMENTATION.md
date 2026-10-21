# Lossless writer implementation status

A standalone GNU Fortran probe used non-advancing formatted reads in 4,096-character chunks and wrote the transferred characters without numeric parsing. It compiled with the installed MSYS2 GNU Fortran 16.1.0. Measured synthetic payloads of 17, 80, 81, 600, and 9,000 characters plus an empty record retained their character content.

The probe normalized line endings to the platform-native output and appended a record terminator when the input's final record had no linefeed. Thus it is **not** the PHITS writer fix required here. Exact CRLF/LF and unterminated-EOF preservation, multi-thread merge order, NCOL=4 schema, NCOL=2 termination, and end-to-end PHITS operation remain untested. No PHITS source was changed, copied into the repository, or compiled. No production patch is offered from this probe.

Probe files and console outputs are local under `D:\temp\gnu_lossless_copy_probe_20261021`; they are not part of the repository.
