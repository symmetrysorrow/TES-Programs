# Run identity and hashes (2026-10-21)

## R100: prior same-run diagnostic dataset, normal CLI conversion

- PHITS input SHA-256: `9624da6b9a2542102b66c251e5d5e232bc2233432357f98dce520563c274c6af`.
- Raw Dumpall: 54,036,774 bytes; SHA-256 `d6f1dad78ed7d768a65a14e6b17642ec7ec14cb14e77d157a85c62ded4a34164`.
- `[t-deposit]`: SHA-256 `91c226eb5d4047b700f9e135ae88d316976dd4165246b8dc596b79550f51d57b`.
- Internal PRE/POST trace: SHA-256 `124d900962e163893010c34c6cc0adfd18e56e881f83ef6fbe1405a80e1aeb6a`.
- Normal CLI HDF5 v2: SHA-256 `6b202db512b8848b4e079588245abff22227f6d340a8d025b2bba9fe574c1eb7`, 11,523,587 bytes.
- This historical run folder records PHITS exit 0. Its exact launch executable hash was not recorded at original execution time; see the candidate diagnostic executable hashes in `PREVIOUS_REPAIR_EVIDENCE.md`.

## R1000: new independent-seed same-run PHITS and normal CLI conversion

- PHITS 3.341 diagnostic executable: `phits_WinGfort_lossless_seq.exe`, SHA-256 `8bb304d858ba1dd4e94436759e8dc38b7f955a80c301e3c20f6c16bc5709c77e`.
- Run input SHA-256: `bdcf9d95c7328735b89a7d3ea0ea4cea40c4a7f1d41eb83a22dc4d27fe5e94ac`; source input before allowed run-local edits: `9624da6b9a2542102b66c251e5d5e232bc2233432357f98dce520563c274c6af`.
- Seed: `1011100101010100011100110010110101110101010111100010110101101011`; maxcas 1,000; maxbch 1.
- Raw Dumpall: 277,977,602 bytes; SHA-256 `9eb5f4abf86b0421ee22e1adbb4f7f78104297c605ceb0f39c8b919c54b4ea`.
- `[t-deposit]`: 11,945 bytes; internal trace: 164,150,585 bytes; exact hashes are in the private run metadata/comparison JSON.
- PHITS exit 0; CPU-reported run completed in 76.74 s.
- Normal CLI HDF5 v2: 58,865,989 bytes; SHA-256 `4dc037236e20bcc9b5a981e84177d1f1aea34cccb94bafe38c69fbe6ba2d7dc4`.

Large raw, trace, table, HDF5, and row-level comparison files remain outside the Git worktree. Their local execution directory was `D:\temp\phits_n1000_diag_20261021\run`.
