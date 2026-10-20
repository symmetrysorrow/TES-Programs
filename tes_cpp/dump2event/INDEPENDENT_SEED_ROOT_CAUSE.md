# Independent-seed histogram difference: root cause

## Classification

**Cause: energy reconstruction difference.** The first raw-only candidate
omitted a PHITS 3.341 NCOL=10 EGS5 electron boundary-crossing scoring branch.
It was not a bin-boundary or floating-point assignment effect.

## Source path

In PHITS 3.341 `src/talls08.f`, `tdepstregEVENT` initializes the energy
variables from the transport record. For a boundary-crossing NCOL=10 record it
retains the EGS5 branch. In `src/egs5edxde.f`, the EGS5 electron stopping
power is evaluated at the incoming energy and computes `ecc=ein-dedx*step`.
The tally contribution in `tdepstregEVENT` is `tlv=ein-ecc`. Therefore this
branch contributes

```text
tlv = dedx(E_in, EGS5 medium) * step
step = Euclidean distance from raw Dumpall (x,y,z) to (xc,yc,zc)
```

The raw Dumpall NCOL=10 records serialize particle ID, ITYP, E, EC, weight,
material/region and both endpoints, enough for this validated homogeneous Pb
segment. In these eight records JCOLL is 0 in Dumpall, although PHITS's
scoring diagnostic identifies the EGS5 state; NCOL and ITYP identify the
validated branch, so the prototype does not require JCOLL=14 here.

## Reproduction

Input Dumpall: existing independent-seed 1,000-history file
`dump2event_physics_validation/overnight_physics_reconstruction/runs/n1000_independent_seed/dumpall.dat`.
The PHITS internal reference is the same-seed diagnostic trace generated from
the matching PHITS 3.341 setup. The raw candidate uses only Dumpall, the local
stopping table and the fixed experimental configuration; it does not read the
internal trace. Compare rows by NOCAS and particle ID, then sum over the
history to reproduce the `[t-deposit]` value.

Initial candidate versus PHITS internal history totals had seven affected
histories, each missing one or two NCOL=10 additions. The eight positive
omissions sum to `0.018778564270170956 MeV`:

| NOCAS | particle ID | PHITS boundary addition (MeV) |
|---:|---:|---:|
| 78 | 2 | 0.00075307754072129196 |
| 78 | 20 | 0.000050214582755354203 |
| 246 | 2 | 0.0034032325973862299 |
| 289 | 2 | 0.0010219964378909152 |
| 365 | 2 | 0.0019244505388539279 |
| 443 | 2 | 0.0029613293557944065 |
| 461 | 2 | 0.0052146536067856353 |
| 871 | 2 | 0.0034496096099831952 |

The detailed record IDs and endpoints are in
`validation/independent_seed_1000_ncol10_boundary_audit.csv`. All eight
reconstructed endpoints matched the internal diagnostic endpoints exactly
(0 cm coordinate difference). Three histories previously crossed a 0.01 MeV
bin edge due to the missing energy: NOCAS 246, 365 and 461. The initial
candidate values were at least `2.3528936882455387e-5 MeV` from the nearest
edge, so these are not floating-point tie cases.

## PHITS histogram rule and final check

PHITS 3.341 `talls08.f` accumulates `trEVENT` for each NOCAS. At history flush
it scores only positive history sums satisfying `heats >= eb(1)` and
`heats < eb(ne+1)`, then assigns the half-open bin `[eb(i), eb(i+1))`. The
configured comparison was `[0,1.5) MeV` at 0.01 MeV per bin. The local tally
weight is `sweight(m)*ratio`; all compared histories have unit scoring weight.

After implementing the NCOL=10 branch, a fresh C++ ledger was compared with
the PHITS internal history trace:

- 1,000 histories, 201,971 tally additions versus 201,971 ledger rows.
- 0 unresolved rows in this selected data set; executable still reports
  `EXPERIMENTAL_PARTIAL` because its physics scope is deliberately limited.
- PHITS internal sum: `35.72084236493933 MeV`.
- Reconstructed sum: `35.720842364939166 MeV`.
- Absolute total difference: `1.6342482922482304e-13 MeV`.
- Maximum per-history absolute difference: `7.815970093361102e-14 MeV`;
  zero histories exceed `1e-12 MeV`.
- Both sides have 58 positive histories in range; histogram L1 difference is
  0 over 150 half-open bins.

The same-seed GNU baseline and diagnostic runs also had byte-identical
Dumpall and `[t-deposit]` output, supporting that the diagnostic logging hook
did not change transport for this run. PHITS's stock distribution output has
no history IDs; the internal trace supplied the IDs for the exact check.

Reproduction scripts are `validation/run_physical_reconstruction.py` and
`validation/compare_internal_history_totals.py`. PHITS diagnostic source,
binary, raw Dumpall and stopping table are local research inputs and are not
included in this branch.
