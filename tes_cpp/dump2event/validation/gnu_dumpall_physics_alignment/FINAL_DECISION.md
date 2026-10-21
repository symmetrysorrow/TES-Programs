# Final decision

**Parser behavior:** the C++ schema parser supports valid wrapped numeric fields and now reports detailed fail-closed context. It cannot restore the GNU ASCII values already truncated by PHITS’s `character(80)` merge buffer. Three GNU files all fail at the first NCOL=4 particle tuple.

**Physics mapping:** this run achieved 0 internal-to-raw addition matches and 0 new physics-value reproductions across the 189,126 internal photon score calls. Existing electron and raw-only positron evidence remains intact but is not evidence from the GNU runs. No new deposit equation or physics branch was implemented.

**Weights/cutoff/annihilation/photon:** aggregate internal diagnostics are retained with provenance; raw matching is unresolved. Per-addition success counts, maximum raw/internal errors and position/region rates cannot be reported because there are no matched pairs.

**Required next action:** regenerate an ASCII Dumpall after correcting PHITS `cpcasdump` to preserve records beyond 80 characters, or use the normal unformatted Dumpall and a schema-aware decoder. That change must be made in a separate diagnostic source/build path, never in `C:\PHITS` or the protected input data. Then rerun mapping before implementing any additional physics.

**Implementation decision:** no new physics behavior is ready for the converter. The parser failure-context improvement is safe and independent. Historical parser and physical-deposit regression results remain documented; this turn did not re-read the unavailable large reference files.
