# PHITS 3.341 same-seed diagnostic run

## Run identity

- PHITS `Version = 3.341` from the run's `phits.out`.
- Run began 2026-10-21 11:57:40.0137401 JST and finished 12:05:25.0059984 JST.
- PHITS output: `END`; batch log: `ncas=10000`, remaining batches 0; exit code file: `0`.
- CPU time reported by PHITS: 7m 37.93s.
- Isolated executable: `phits_WinGfort.exe` from a temporary diagnostic source copy. The installed PHITS tree was not built into or modified by this diagnostic change.
- Input path recorded by the run: `tes_cpp\dump2event\validation\positron_10000_closure\trace_10000_diagnostic_template.inp`, executed from an isolated temporary output directory. The input hash is `7401c49903c41317825af7eed2a8c239033ea7505514cecebe0f1524459f007e`.
- Input keeps `icntl=0`, `maxcas=10000`, `maxbch=1`, source seed, 1.173/1.332 MeV spectrum, `negs=1`, `e-mode=2`, `emin(12/13/14)=0.001`, model/material/geometry, and Pb region 1000 `[t-deposit]`; only `file(15)` redirects the newly generated Dumpall to the isolated run directory.
- The new-run Dumpall was 2,736,869,598 bytes and SHA-256 `c11c51299c6c77a7e8a531a640ed72b31c8455fb2d0a17f2b03d3716a2815060`. It contains 10,000 histories and ends with NCOL=2.
- The baseline and diagnostic `[t-deposit]` output SHA-256 is identical: `0d87d8b7041c54cfd10e610ac0653ff1eda1844a31ca0691126ce82ff0c5f110`.
- Diagnostic `tally_internal_trace.log` contains 1,954,742 `PRE` score-addition rows; output size 1,643,938,022 bytes. It remains local/private.

## Command recording limitation

The process metadata preserves the executable path, input path, working directory, process exit code, run start/finish time, and complete PHITS/batch logs, but the parent shell did not persist the literal `CreateProcess` command line. It would be inaccurate to present a reconstructed command as the exact executed command. Re-run command form for this isolated executable is:

```powershell
& '.\phits_WinGfort.exe' '.\trace.inp'
```

with the run directory as the working directory. This command form is supplied for reproduction and is not claimed as a byte-for-byte captured invocation.

## Output directory and safety

All outputs (PHITS log, Dumpall, internal trace and t-deposit history) were written under a dedicated temporary diagnostic run directory. Original 10k Dumpall size, modification time and SHA were unchanged. The 116.8 GB original Dumpall, original PHITS input, PR #6 and installed PHITS tree were not written.
