# Previous GNU Dumpall repair evidence search

Date: 2026-10-21. Search roots included `D:\temp`, the research worktree, and `C:\Users\Nao\Documents` / `Downloads`; no prior PHITS source patch, repaired-source tree, or repaired PHITS executable was found in those locations. Existing records are reports and outputs only. The previous reports identify one 100-history repaired Dumpall (`n100_lossless_v2`, 54,036,774 bytes, SHA-256 `d6f1dad78ed7d768a65a14e6b17642ec7ec14cb14e77d157a85c62ded4a34164`) but do not identify the patch source or corresponding binary hash. That artifact itself was not found at its historical path in this environment.

A standalone GNU Fortran 16.1.0 line-copy probe was compiled under `D:\temp\gnu_lossless_copy_probe_20261021`. Content checks passed for 17, 80, 81, 600, and 9,000 character records and an empty record. The probe emits platform-native formatted record separators and adds a line ending to an unterminated last line. It is a primitive test only, not a PHITS merge repair; it fails the requested exact EOF/newline preservation test and does not prove thread ordering or NCOL structure.

**Status:** prior PHITS repair implementation and binary identity not recovered; no repair claimed.
