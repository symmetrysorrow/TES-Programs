# Scoring completeness versus full absorption

The HDF5 event `complete` field means only that all score-candidate rows encountered for that event were reconstructed (no UNRESOLVED rows). It does not claim that all transport physics is represented in Dumpall, that the energy balance closes, or that the event was fully absorbed.

The converter additionally requires a structurally complete Dumpall ending in NCOL=2. Provenance labels transport physics completeness `undetermined`, energy balance `not assessed`, and all absorption `undetermined`. Leakage (including NCOL=12) remains visible in tracks and is not erased by per-event renormalization. TES block assignment and pulse generation remain downstream event2pulse responsibilities.
