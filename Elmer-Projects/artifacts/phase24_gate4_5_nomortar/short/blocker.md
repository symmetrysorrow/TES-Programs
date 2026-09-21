# Gate 4/5 blocker: no-mortar on production-v2 nonconformal mesh

The no-mortar steady current reference is numerically close to COMSOL, but it
does not provide a valid transient initial field on `mesh_singlepixel_prod_v2`.
The mesh is nonconformal at the interface; disabling mortar leaves the
interface coupling invalid for the transient restart path.

Evidence from the no-mortar MUMPS Gate 4 short-window attempt:

- Case: `case_phase24_gate4_5_nomortar_short_mumps`
- Restart: `case_tes_steady_prod_v2_nomortar.result`
- Interface: `Apply Mortar BCs = False`
- `TESInnerCircuit` temperature after restart: about `-1.36e5 K`
- Heat `Result Norm`: about `22625.86` (the valid mortar transient is about `0.17`)
- The case was stopped before waveform comparison; Gate 4 and Gate 5 are not
  passed by this bundle.

Conclusion: the arithmetic Gate 3 current check is insufficient for this
nonconformal no-mortar route. A valid Gate 4/5 campaign requires either a
conformal no-mortar mesh or restoring mortar for the production-v2 mesh.
