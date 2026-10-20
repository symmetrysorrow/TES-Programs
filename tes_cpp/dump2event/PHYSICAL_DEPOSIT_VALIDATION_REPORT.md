# PHITS physical-deposit reconstruction validation

## Scope and result

This is a narrow research prototype for PHITS 3.341, Pb region 1000 (material
ID 2 / EGS5 medium 2), unit weights, 1 keV electron cutoff, and the matching
PHITS-derived stopping-power coefficients. It writes a row ledger only. It is
not a complete event converter and must not be treated as general PHITS
physics support.

The main independent-seed issue was traced to omitted NCOL=10 EGS5 electron
boundary loss. Adding the source-derived term makes every history total match
the same-seed PHITS internal tally trace within floating-point accumulation
error. Quantitative results are recorded in
`INDEPENDENT_SEED_ROOT_CAUSE.md` and machine-readable files under
`validation/`.

## Retained baseline checks

The existing 100-history raw-only validation reports 38,077/38,077 additions
matched by ID, NOCAS and region, with maximum absolute energy error
`5.7976626971489864e-15 MeV`, no positive coordinate mismatch, and all 38,077
rows accounted for (37,190 continuous terms, 757 other positive terms, and
130 exact-zero rows). The 10,000-history structural parser read 1,983,461
records and reached NCOL=2. These are retained validation results, distinct
from the independent-seed comparison.

## Independent-seed exact comparison

| Metric | Result |
|---|---:|
| Histories | 1,000 |
| PHITS internal tally additions | 201,971 |
| Raw Dumpall ledger rows | 201,971 |
| Unresolved rows in this sample | 0 |
| PHITS total | 35.72084236493933 MeV |
| Reconstructed total | 35.720842364939166 MeV |
| Absolute total difference | 1.6342482922482304e-13 MeV |
| Maximum absolute history difference | 7.815970093361102e-14 MeV |
| Histories with error above 1e-12 MeV | 0 |
| Positive histories in `[0,1.5)` | 58 on both sides |
| Histogram | 150 bins, 0.01 MeV, half-open; L1 difference 0 |

The comparison is keyed by PHITS internal NOCAS and compares sums before
binning. It does not infer history identities from the ordinary PHITS
histogram. Per-history values and bin edges are in
`validation/independent_seed_1000_history_compare.csv` (generated during
local validation; not committed if it is too large for this branch).

## Current implementation rules

- NCOL=13/14, JCOLL=14, ITYP=12: `dedx(E) * endpoint chord`.
- NCOL=10, ITYP=12: `dedx(E) * endpoint chord`, matching the EGS5 boundary
  scoring branch.
- NCOL=11, ITYP=12: residual electron energy plus `dedx(1 keV) * final
  chord`, for this specific validated cutoff setup.
- Dead electron secondary above the conservative `1e-6 MeV` score guard:
  secondary kinetic energy at its serialized point, for unit secondary weight.
- Exact `E == EC` in supported ordinary score records: an exact zero with no
  point claim.

No negative clipping, global renormalization, or arbitrary energy transfer is
performed. Unsupported rows stay unresolved. The low-energy guard is
conservative because `esmin` is not recorded in Dumpall. This prototype does
not yet prove segment splitting at region boundaries or voxel-level energy
placement.

## Build and tests run

On 2026-10-20, MSVC 19.51 / Visual Studio 18 built `dump2event`,
`dumpall_record_inspect`, `dump2event_physical_reconstruct`, and
`test_physical_deposit_engine`. CTest passed 1/1. Python unittest passed
23/23 with the built converter selected. The default converter processed the
original 1,000-history sample with `--summary-only --save-all` and reported
1,000 histories, exit code 0. No Linux/WSL build was run.

## Decision

Research use is reasonable only for the explicitly stated PHITS 3.341 / Pb /
unit-weight validation slice and with all unresolved rows reviewed. It is not
ready for PR #6: no general scoring branch coverage, positron rest-mass
accounting, arbitrary material/weight handling, voxel allocation, or full
reconstruction of all PHITS tally modes has been validated.
