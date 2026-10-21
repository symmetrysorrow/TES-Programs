# PHITS run inputs, executable and outputs

PHITS version was read from each `phits.out` as `Version = 3.341`. The exact-sequence serial executable was `phits_WinGfort_lossless_seq.exe`, SHA-256 `8bb304d858ba1dd4e94436759e8dc38b7f955a80c301e3c20f6c16bc5709c77e` (44,143,219 bytes). It ran in isolated per-run output directories; PHITS installation, repository input, and source template were not edited. The run wrapper changed only `maxcas` to the requested number and retained `maxbch=1`, seed, geometry, material, source, cutoffs, and tally cards.

## Completed 100-history runs

| Run | Start–finish (JST) | Derived input SHA-256 | Dumpall bytes / SHA-256 | `[t-deposit]` SHA-256 | PHITS exit |
|---|---|---|---|---|---:|
| 0.511 MeV | 2026-10-21 18:43:54–18:45:08 | `a813f4f01f84ce78dd72491e39b0351fcdac792011af4974307cd74e5804a018` | 354,930,838 / `1f3948816f016daa0fb6de004e7def521ad8480e0282b38037774d9a7c15fef8` | `b654ec761b89daffab4145bc0c20d89bc911f6777823e5bf8edcb888f887e854` | 0 |
| 2 MeV | 2026-10-21 18:45:08–18:47:49 | `5a6bb474879c2cb0a84b9811f4f744206654ced823dada1f335f5d32237bea5a` | 777,546,730 / `41da591838063bd9789fcc9994060aa3c470d00939b116b79aaf5ecbe171d704` | `7eb262aaed6396a4277897834264348d83f46d729c32fb5069bd2046a0bccd6d` | 0 |
| 0.511 MeV / 1,000 | 2026-10-21 18:57:59–19:14:28 | `f84fe3dd252807064dea6b8ad255e472260de6ff14ef2992c1ee5e98694c4466` | 3,594,947,832 / `4d06d8311232b1ba8c39e32e5341460d9fdaa187d1f6c4c7fe2e760e29df36f8` | `9000d9dcd65e6d7fce152493fd4d231554a16798b1e52bc94d74035a822e5c47` | 0 |
| 2 MeV / 1,000 | 2026-10-21 18:58:13–19:33:58 | `1d6bc349f47533bf6f4eca1835313977c858050b233466e0e3d32add74f86b3c` | 8,457,662,268 / `f84d0dd4624549bfcc1b54e705109aecc93c3b9a7ee345089d16d49baabc59dd` | `4823687fcc3b9927f02e9a4912a693522f0cfba163016dca68fe690c5ccd4b0e` | 0 |

Each PHITS stdout and `batch.out` reported no remaining batches at termination. The full run-local trace, sequence sidecars, input, output, and file hashes are kept outside the Git worktree. `[t-deposit]` and Dumpall comparisons always use the same run directory.

The diagnostic command for both 100-history runs was the repository's isolated PowerShell runner with `-InputFile <mono-gamma template> -RunDirectory <new run dir> -MaxCas 100 -Executable <lossless serial PHITS executable>`. It copied the input into the new directory, changed only `maxcas`, set `PHITSPATH=C:\PHITS` and `OMP_NUM_THREADS=1` for that child, fed `file=trace.inp` on stdin, and captured stdout/stderr and the exit code. The corresponding 1,000-history commands used the same runner and inputs with `-MaxCas 1000`.

The normal conversion command was `dump2event.exe <same-run Dumpall> <new event.h5> --physics pb-egs5 --hdf5-version 2 --physics-config <run contract JSON> --stopping-table <Pb EGS5 coefficient CSV>`. `--normalize-deposits` was omitted. For the default-mode compatibility smoke test, the command was `dump2event.exe <511-keV raw Dumpall> <new event.json> --input-energy 0.511`; it exited 0 and wrote 60 retained events.

## 1,000-history runs

Both 1,000-history diagnostic jobs used the same mono-gamma definitions and 64-bit seed as their corresponding 100-history runs; only `maxcas` was set to 1,000. The 0.511 MeV job exited 0 and produced a 3,594,947,832-byte losslessly written Dumpall, SHA-256 `4d06d8311232b1ba8c39e32e5341460d9fdaa187d1f6c4c7fe2e760e29df36f8`. Its CLI reported 2,541,373 raw records, 1,000 histories, terminal NCOL=2, and 2,541,564 HDF5 deposits with 0 negative, 0 unresolved, and 0 unsupported rows. The exact same-run join completed for all 2,541,564 additions; the maximum energy error was 4.119512890532473e-14 MeV, the maximum event-sum error was 1.1268763699945339e-13 MeV, and region/identity/position mismatch counts were zero. Its 150-bin tally comparison was exact.

The 2 MeV / 1,000 PHITS job finished at 2026-10-21 19:33:58 JST with exit code 0 after 35m39.27s and `remaining batches=0`. Its PHITS output SHA-256 is `6ac234ac2834306aa3a0c8bf2c5d0daa70921909bf81f61ad964111854cfe0ba`. `stderr.log` contains the runtime notice `IEEE_INVALID_FLAG`; it did not change the process exit status or stop batch completion, and it is retained as a diagnostic caveat. Its Dumpall parser reported 5,978,725 records, 1,000 histories and terminal NCOL=2. The normal CLI withheld HDF5 because one dead secondary is below the validated score guard. Independent raw-ledger comparisons nevertheless matched all 5,978,912 reconstructed increments and all 150 histogram bins to the same-run PHITS tally; one candidate remains explicitly unresolved. Details are in `UNRESOLVED_PHYSICS.md` and `INTERNAL_ADDITION_RECONCILIATION.md`.

## 10,000-history CLI attempt

The existing complete Co-60 10,000-history input generated `dumpall.dat` (2,736,869,598 bytes, SHA-256 `c11c51299c6c77a7e8a531a640ed72b31c8455fb2d0a17f2b03d3716a2815060`) using PHITS 3.341 and exit code 0. A new normal CLI attempt failed at the first NCOL=4 header because the legacy logical row was malformed for the strict parser. No HDF5 was published. The original Dumpall and PHITS input remain unchanged; see `TEN_THOUSAND_HISTORY_CLI_REPORT.md`.
