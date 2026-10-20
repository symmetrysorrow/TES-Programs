# GNU Fortran Dumpall NCOL=4 root cause

Date: 2026-10-20 (Asia/Tokyo). PHITS source inspected: `C:\PHITS\src\ovly15.f` from the installed PHITS 3.341 tree. Inputs were read only.

## Finding

The NCOL=4 record is emitted by `dumpdat` as separate list-directed writes. The history header is followed by a particle tuple of eight values in this order: `no, idmn(mat), ityp, ktyp, jtyp, mtyp, rtyp, oldwt`. This is explicit at `ovly15.f` around lines 1064–1083. Energy, weight, position, regions and optional `qs` are later independent write statements.

In the observed GNU-built diagnostic Dumpall, all copied ASCII records are exactly 80 characters. The first NCOL=4 particle tuple ends at character 80 as `... 0.000`; the next physical line is the four-value region tuple. Therefore the 80-character intermediary has truncated a logical write, leaving fewer than the eight required values. This is data loss, not merely a token wrapping convention.

The source path is `dumpdat` → `cpcasdump` when a new history begins. `cpcasdump` declares `character(80) dummy` (line 1328) and, for `idumpall=-1`, copies records using `read(io2,'(a)') dummy` followed by `write(io,'(a)') dummy` (lines 1486–1490). A source record longer than 80 characters cannot survive this copy. The 80-character file records and clipped numeric field directly match that code path.

Observed first five physical lines are identical in structure in all three GNU diagnostic files:

| File | Size (bytes) | NCOL=4 particle line | Following line | Parser result |
|---|---:|---|---|---|
| `D:\temp\weighted_cutoff_positron_photon_20261020\diag_gamma2_20\gamma_2mev_thickPb_dumpall.dat` | 172,131,654 | 80 chars; 7 tokens, last token `0.000` | 4 region values | 11 combined values, expected 8 |
| `...\diag_gamma511_thick_20\gamma_511kev_thickPb_dumpall.dat` | 73,544,120 | same | same | same structural failure |
| `...\diag_gamma511_thin_20\gamma_511kev_thinPb_vacuum_dumpall.dat` | 17,330,578 | same | same | same structural failure |

The underlying raw temporary per-thread file was not retained, so its full unclipped bytes could not be compared. The exact values lost from the persisted files are not inferable. No values were imputed.

## Outcome

The parser now accepts legitimate list-directed values split over physical lines when the schema count is satisfied, and emits source path, physical line, record start, NCOL and expected field context for the observed malformed NCOL=4 tuple. It deliberately rejects these GNU files. Regeneration through a corrected PHITS ASCII merge routine or an untruncated binary Dumpall is needed for record-level physics mapping. No PHITS source or installation file was changed.
