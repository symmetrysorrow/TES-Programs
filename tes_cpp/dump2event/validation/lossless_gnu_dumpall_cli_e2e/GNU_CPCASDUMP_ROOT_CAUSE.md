# GNU `cpcasdump` truncation: root cause

Date: 2026-10-21 (Asia/Tokyo). PHITS installation source was inspected read-only at `C:\PHITS\src\ovly15.f`.

`dumpdat` writes the NCOL=4 particle tuple with eight values. In `cpcasdump`, the GNU ASCII merge path declares `character(80) dummy` at line 1328 and copies each formatted record with `read(io2,'(a)') dummy` / `write(io,'(a)') dummy` at lines 1486–1490. The observed GNU diagnostic outputs contain exactly 80-character lines; the NCOL=4 tuple ends with a clipped `0.000` fragment and the next region tuple yields 11 values where the parser schema expects 8. The lost digits cannot be recovered from those files.

The prior finding is confirmed by installed source and raw line shape. No attempt was made to infer missing fields. The `C:\PHITS` source and binaries remain unchanged.

**Status:** root cause identified; production PHITS writer not yet repaired or rebuilt in this task.
