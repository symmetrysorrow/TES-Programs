# Commands and measured checks

Date: 2026-10-21 (Asia/Tokyo). This records commands executed for this task; no new PHITS transport run was executed.

| Command / action | Result |
|---|---|
| `git -C D:\github\TES-Programs-physdeposit-20261020 rev-parse HEAD` | `b19d640e9fb6b7fd0e6c2cf4dda3153af57670ca` |
| `git -C ... branch --show-current` | `feat/dump2event-physical-validation-20261020` |
| `git -C ... rev-parse origin/feat/dump2event-physical-validation-20261020` | same SHA as local HEAD |
| `cmake --build build/physdeposit-vs --config Release -j 4` | exit 0 |
| `ctest --test-dir build/physdeposit-vs -C Release --output-on-failure` | exit 0, 1/1 passed |
| `python -m pytest tes_cpp -q` | exit 0, 44 passed, 10 skipped |
| `run_copy_probe.ps1` with MSYS2 GNU Fortran 16.1.0 on PATH | exit 0; payload preserved in all seven cases; final LF was added for the no-final-LF case |

The script is in `scripts/run_copy_probe.ps1`; it builds and runs a generic Fortran record-copy probe. It does not invoke PHITS. The GNU test observed input/output byte counts: short 18/19, exact80 81/82, long81 82/83, hundreds 601/602, multi_k 9001/9002, empty 1/2, no_final_lf 81/83. On Windows, output record terminators are CRLF. Character payloads matched; exact line-ending preservation did not.

Existing historical PHITS evidence remains as described in the individual reports. The HATCH table, repaired PHITS executable, lossless Dumpall, HDF5 v2 result, and a same-run raw/internal join were not obtained in this task.
