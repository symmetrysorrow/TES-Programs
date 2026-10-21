# Run identity and SHA audit

Audit date: 2026-10-21. No new PHITS run was executed in this task.

| Historical artifact | Size / SHA-256 | Audit status |
|---|---|---|
| Earlier 10,000-history raw | 2,783,824,000 / `d254b3afd13230be84ad8f733c8aaa89811724f832a506a0d5b421d051c06b1c` | Parser accepted 1,983,461 records and NCOL=2 terminal; not same run as Oct 21 diagnostic trace. |
| Oct 21 baseline raw | 2,736,869,598 / `c11c51299c6c77a7e8a531a640ed72b31c8455fb2d0a17f2b03d3716a2815060` | GNU ASCII output; parser rejects first NCOL=4 particle tuple because values were truncated. |
| Oct 21 diagnostic raw | 2,736,869,598 / same `c11c51299c6c77a8e8a531a640ed72b31c8455fb2d0a17f2b03d3716a2815060` | Byte-identical to baseline; malformed for record-level mapping. |
| Historical 100-history lossless-labeled raw | 54,036,774 / `d6f1dad78ed7d768a65a14e6b17642ec7ec14cb14e77d157a85c62ded4a34164` | Reported previously; raw and matching writer executable/source were not located for this task. |

The exact all-artifact same-run closure requested by this task was not achieved. None of the historical internal trace or tally files is joined to the malformed October GNU raw as if it were valid.
