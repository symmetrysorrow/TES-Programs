# Negative and duplicate audit

No negative, NaN, or infinite deposit was written in any output. The v2 writer rejects negative and non-finite energies; all three stored HDF5 files were reread through structural/offset validation. Addition IDs are checked for duplicates both within an event and across prior events; event and particle IDs are also checked. The 100-history HDF5 was additionally keyed against the internal addition oracle with no duplicate, missing, or extra IDs.

| Run | Negative HDF deposits | Unresolved deposits | Duplicate-ID/offset integrity |
|---|---:|---:|---|
| 100 | 0 | 0 | pass; 38,077/38,077 IDs match oracle |
| 1,000 | 0 | 0 | pass; HDF5 structure and offsets valid |
| 10,000 | 0 | 342 | pass; unresolved rows have invalid energy, not zero |

No negative-value clipping, history renormalization, or fallback to the parent-minus-secondary formula was used.
