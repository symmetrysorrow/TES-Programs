# Pb HATCH table availability, rights review, and installation behavior

Date: 2026-10-22 (Asia/Tokyo)

## Verified source bytes

The original file exists in the local, untracked validation directory:

`dump2event_physics_validation/egs5_stopping_power_reconstruction/runs/n100_stopping_probe_v2/egs5_stopping_power_coefficients.csv` relative to the main checkout

Measured values:

- Size: 130,754 bytes
- SHA-256: `20fd86eb0b57483bd19f20a8388d2c1ecb0431a02fa99f0933e2c3c4748514e9`
- CSV header: `medium,medname,lelke,eke0,eke1,ededx0,ededx1,pdedx0,pdedx1,phits_density,egs_density`
- Data rows: 600 total; medium 2/PB: 150
- Medium 2 indices: unique `lelke` 1 through 150
- PHITS and EGS5 density values: 9.78 g/cm3

The original file was read only. No coefficient values were changed, reconstructed, or substituted. The similarly named `n100_stopping_probe` copy was not selected.

## Public redistribution review

**Decision: redistribution permission is not established; the CSV is not included in Git.**

The PHITS official download/licensing page states that PHITS is supplied under an individual-use license and each user must obtain a license. It does not grant a public redistribution right for HATCH-generated output data. The installed `C:\PHITS\readme-eng.txt` points to license agreements in the distribution's `document` directory. A filename inventory of its 50 files found no PHITS license/agreement document; the `C:\PHITS\document\license` subdirectory contains third-party tool licenses. No PHITS output-data redistribution clause was found locally.

The EGS5 manual's SLAC notice says the manual and included software are for non-commercial use and that notices must remain attached to copies/derivatives of the manual/software and documents/software containing such material. That notice does not specifically classify or grant redistribution rights for HATCH-generated coefficient CSV output. Whether this table is covered as a derivative/output is therefore unresolved. We did not infer permission merely because it was produced by a local run.

Sources:

- PHITS official Japanese licensing/download procedure: https://phits.jaea.go.jp/getj.html (individual license condition)
- EGS5 official KEK-hosted manual, SLAC Software Notices, printed page i / PDF page 1: https://rcwww.kek.jp/research/egs/egs5_manual/slac730-150228.pdf
- Local PHITS package note: `C:\PHITS\readme-eng.txt`, lines describing `document` as containing updated-feature documents and license agreements

A written clarification from the relevant PHITS/EGS5 rightsholders would be needed before the CSV can be added to the public repository.

## Code and local use

The resource resolver now enforces the historical SHA-256 for all routes, including an explicit `--stopping-table` path. It rejects arbitrary, altered, or wrong-table CSVs. Search order without an explicit path is:

1. `TES_CPP_PB_HATCH_TABLE` (private per-user installation)
2. CSV beside the Python driver (setuptools/Python package install)
3. source-tree `tes_cpp/dump2event/data/`
4. `<CMake install prefix>/bin/`

The installed CMake layout previously placed the driver under `share/tes_cpp_dump2event` while binaries/resources belong in `bin`; the resolver now checks the prefix `bin` directory. `setup.py` copies the driver beside the executable and can copy a table only if one has been legitimately installed into the source resource directory.

PoST's interactive flow checks the same user environment setting, validates its checksum before skipping the table picker, and rejects non-Pb material. It does not offer Bi as a usable path. HDF5 v2 provenance continues to store the SHA-256 of the exact table passed to the converter.

Example once-per-user setup:

```powershell
[Environment]::SetEnvironmentVariable(
  "TES_CPP_PB_HATCH_TABLE",
  "D:\private\egs5_stopping_power_coefficients.csv",
  "User")
```

Then the requested normal command is:

```powershell
dump2event input.dumpall output.h5 --physics-config pb_config.json
```

No coefficient CSV was copied into this worktree or staged for Git.

## Existing data and regression checks

Read-only inspection of the existing the existing `gamma_2mev_1000.h5` in the local `two_mev_1000` validation directory found HDF5 v2, PHITS 3.341, algorithm `pb-egs5-exp-v3-score-gate`, and the same stopping-table SHA-256. The file was not regenerated or modified.

The archived report for this file records 5,978,912/5,978,912 internal additions matched, 0 addition/event/region/coordinate mismatches, and maximum energy difference `4.680023729664029e-14 MeV`. This resource task changed only resolution and UI behavior, not physics equations; that recorded physical validation remains the prior evidence.

New software checks cover checksum rejection for explicit and environment table paths, missing/corrupt resources, source/package path discovery, and Pb-only UI flow. Results:

- MSVC Release CMake build: success with MSVC 19.51.36260.0
- CTest (`tes_cpp/build/dump2event`): 1/1 passed
- pytest (`tes_cpp/dump2event/tests`): 77 passed, 7 skipped
- `pip install . --no-deps --target <temporary-install-directory>`: success; generated platform wheel `tes_programs-0.1.0-cp314-cp314-win_amd64.whl` (281,606 bytes)
- Installed package smoke check: CLI, both physical helper binaries, and the Python driver were present; with `TES_CPP_PB_HATCH_TABLE` set to the exact private CSV, the installed driver resolved and SHA-verified it.
- The first packaging attempt used `--no-build-isolation` and failed because the global Python 3.14 environment has no `setuptools.build_meta`. The isolated pip build then succeeded without installing build dependencies globally.
- Existing `gamma_2mev_1000.h5` was read without rewriting: HDF5 v2, PHITS 3.341, and table SHA-256 match. The stored audit report records 5,978,912/5,978,912 additions matched (max difference `4.680023729664029e-14 MeV`, zero ID/region/coordinate mismatches).

## Git state and merge gate

Work is isolated in `D:\github\TES-Programs-pb-hatch-20261022`, branch `feat/pb-hatch-resource-20261022`, based on `origin/main` commit `4582f6441ed688ddacd5b96b6dc483ce7af6a7a3`. Main checkout and the existing research worktree were left untouched. PR #6 was not edited. The verified CSV, PHITS distribution/source, Dumpall, PHITS inputs, diagnostic traces, and HDF5 files are not included.

Because the CSV is absent, local installation still requires the existing private file to be located via `TES_CPP_PB_HATCH_TABLE` (or explicit path); true out-of-box conversion without that per-user resource is not possible until redistribution is authorized.
