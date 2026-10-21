# Existing electron regression baseline

The previously verified evidence remains documented in `lossless_gnu_dumpall_cli_e2e/`:

| Regression | Internal additions | HDF5 deposits | Max absolute addition error | Position/region/identity mismatch | Histogram L1 |
|---|---:|---:|---:|---:|---:|
| R100 | 38,077 | 38,077 | 5.7976626971489864e-15 MeV | 0 | 0 |
| R1000 independent seed | 195,185 | 195,185 | 5.830839283627043e-15 MeV | 0 | 0 |

The current MSVC Release build and C++ engine test still pass. Python tests include the provenance additions. A fresh R1000 conversion through the newly rebuilt provenance-aware CLI is running against the same 1,000-history raw file; it must finish before the revision-2 regression is marked complete. The original R100 raw file is not present in this isolated worktree's external run area, so its published prior result is preserved as a regression reference rather than relabeled as a new run.
