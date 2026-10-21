# Addition count closure

The prior 10,000-history mapping reports 1,954,742 internal calls and 1,954,742 candidate ledger rows, but its raw/internal run lineage is under audit (see `VALIDATION_PROVENANCE_AUDIT.md`). Counts below describe the stored candidate ledger classifications, not newly established independent-seed validation.

| Exclusive candidate class | Rows |
|---|---:|
| EGS5 electron continuous | 1,901,831 |
| Electron boundary | 83 |
| Electron cutoff | 38,266 |
| Dead electron secondary | 1,310 |
| EGS5 positron continuous, including zero-valued zero-chord candidates | 1,795 |
| Dead positron secondary | 9 |
| Exact zero, remaining classes | 11,448 |
| **Total** | **1,954,742** |

The old summary added 1,802 repaired rows to the old categories, yielding 1,954,740. The difference is explained at the ledger-classification level by two e+ continuous candidate rows that were already counted among the earlier exact-zero rows and are reclassified in the new exclusive breakdown:

| Candidate addition ID | NOCAS | particle | source line | NCOL/JCOLL/ITYP | E=EC (MeV) | raw chord | candidate reconstruction |
|---:|---:|---:|---:|---|---:|---:|---:|
| 96302 | 507 | 2 | 1,644,802 | 13/2/13 | 0.00097296705103067094 | 0 cm | 0 MeV |
| 1,205,287 | 6121 | 3 | 20,578,897 | 13/2/13 | 0.00096347151926767793 | 0 cm | 0 MeV |

Both have exactly coincident raw endpoint coordinates; the candidate PDedx-times-chord expression is zero. The existing internal mapping has a same-source-line dead-secondary entry with the positron kinetic energy, but it does not expose a separately auditable zero-tlv primary row for these candidate IDs. Accordingly the **ledger class and candidate value are identified, while direct internal-tlv-to-these-two-candidate-ID attribution remains unverified**. Do not interpret the row-count arithmetic alone as full addition identity.

The new exclusive count closure is checked by `test_addition_count_closure.py`. Once a same-run repaired diagnostic raw file is available, rerun the all-additions comparison and replace the provisional lineage statement with direct ID-level evidence.
