# Unsupported and unresolved physical cases

The experimental executable is intentionally scoped to a specific PHITS
3.341 Pb validation. The following are not validated as physical deposit
reconstruction:

- Any PHITS release other than 3.341, region/material/EGS5 medium other than
  Pb region 1000 / material 2 / medium 2, or non-unit parent/secondary weight.
- Photon collisions and photon transport losses; positron slowing, stopping
  and annihilation, including rest-mass energy and annihilation photon escape.
- Nuclear, neutron, heavy charged particle, kerma, dose, and other tally
  scoring branches.
- Dead charged secondaries other than the validated electron case, or
  secondary kinetic energy `<= 1e-6 MeV`. The tally threshold depends on
  `esmin`, which is not serialized in Dumpall; such candidates remain
  unresolved instead of being silently counted or dropped.
- NCOL=11 stopping records for non-electrons. They are emitted as unresolved
  even when raw E and EC compare equal because this prototype does not infer
  the scoring branch from equality alone.
- Transport crossing multiple regions/materials within a step, exact
  intra-segment deposition density, boundary splitting, and spatial voxel
  assignment. Endpoint chords are not a validated microscopic point
  deposition.
- General source weights, variance reduction, tally ratios, time/species
  options, and histories where tally normalization differs from unit weight.
- Full event JSON/HDF5 output, production CLI integration, and automatic
  representation of unresolved energy in the existing converter.

The program labels every run `EXPERIMENTAL_PARTIAL`. Unsupported candidate
rows have unresolved status/reasons; no clipping, global normalization,
arbitrary redistribution, or replacement by zero is used. A zero unresolved
count in the 1,000-history sample means only that this particular sample did
not exercise an unsupported row in the selected tally scope.

## Weight accounting confirmed from PHITS 3.341 source

In `src/talls08.f:tdepstregEVENT`, the validated EGS5 `itdfn=0` deposit path
forms `tlv=ein-ecc` and adds that energy to `trEVENT` without multiplying the
sum by a particle weight. At each scored contribution PHITS stores
`sweight(m)=oldwt`; if a later scored contribution in the same history has a
different `oldwt`, it emits the source weight-change warning. At the NCOL=0/4
history flush, the tally score is `sweight(m)*ratio`, with `ratio` determined
by the `part=all` / tally mode branch. Therefore a changing-weight history
cannot be reconstructed by merely multiplying each raw energy term by that
term's `OLDWT`, and a unit-weight result does not validate splitting,
roulette, or importance handling. The 1,000-history comparison described in
the status report exercised unit weights only.

## PR #6 integration blockers

Before production integration, the implementation needs validated
positron-annihilation/rest-mass treatment, photon and other particle score
branches, non-unit weights, exact cutoff semantics, region-boundary segment
splitting, voxel-level comparisons, and a stable runtime provenance/config
interface. It also needs a complete event-output design that preserves every
unresolved contribution without implying completeness. The independent-seed
validation alone does not close these gates.
