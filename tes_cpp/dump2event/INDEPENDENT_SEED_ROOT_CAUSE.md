# Independent-seed histogram difference: root cause

## Classification

**Cause: energy reconstruction difference.** The first raw-only candidate
omitted a PHITS 3.341 NCOL=10 EGS5 electron boundary-crossing scoring branch.
It was not a bin-boundary or floating-point assignment effect.

The `2.22e-15 MeV` code-to-code maximum quoted in the task is not the relevant
PHITS discrepancy: agreement between two raw-only implementations cannot
detect an energy term that both omit. The committed PHITS-internal comparison
records the corrected C++ ledger against the same-seed internal tally, with a
maximum history difference of `7.815970093361102e-14 MeV`. The pre-correction
comparison had eight omitted positive additions and a histogram L1 count
difference of 6. The exact pre-correction raw histogram is not committed; the
three bin crossings below are reconstructed from the archived corrected
history totals minus the individually archived omitted additions.

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

PHITS 3.341 `src/talls08.f:tdepstregEVENT` accumulates each raw `tlv` into
`trEVENT` in Dumpall record order (the `trEVENT(...)=trEVENT(...)+tlv`
statements in the NCOL>=9 scoring branches). On NCOL=0/4 history flush it
selects `heats=trEVENT` (or the configured probability tally path), accepts
only `heats > 0`, `heats >= eb(1)`, and `heats < eb(ne+1)`, then assigns the
half-open interval by direct comparisons `heats >= eb(ie)` and
`heats < eb(ie+1)`. There is no rounding-to-nearest operation before bin
selection. The output records serialize 150 contiguous linear edges from
0.00 to 1.50 MeV (0.01 MeV nominal width); these serialized edges, rather than
an inferred center convention, were used for the cross-check. Underflow,
overflow, zero, and negative history totals do not enter a bin. The histogram
score is `sweight(m)*ratio`; the compared histories are unit weight and use
`part=all`, so the effective weight is 1.

The scorer sets `sweight(m)=oldwt` whenever a deposit addition is made and
prints its weight-change warning if `oldwt` differs from the preceding scored
addition. The physical `tlv` accumulation is separate from this histogram
weight. Thus this exact unit-weight case does not validate mid-history
splitting/roulette or a changing weight. Repeated floating-point additions
follow Dumpall record order; the C++ ledger groups by event and can differ by
roundoff from PHITS' Fortran accumulation order.

Before the NCOL=10 correction, the seven affected histories had eight missing
boundary terms totalling `0.018778564270171039 MeV`. Using the PHITS-serialized
bin edges and the archived current history totals minus those omitted terms,
the old and internal bins differ for exactly three histories: NOCAS 246
(0.8399032761349978 -> 0.8433065087323753 MeV, bin 83 -> 84), NOCAS 365
(0.1091828168095261 -> 0.1111072673483944 MeV, bin 10 -> 11), and NOCAS 461
(0.05665439916464846 -> 0.06186905277144113 MeV, bin 5 -> 6). Each crossing
moves one history from a bin to its neighbor, hence L1 `3*2 = 6`. NOCAS 78,
289, 443, and 871 changed energy but stayed in the same bin. The nearest
pre-correction history to a boundary was still `2.352893688215181e-5 MeV`
away; the nearest corrected positive history was `4.738516994173825e-5 MeV`
away. The observed L1=6 is therefore accounted for by omitted energy, not
floating-point tie handling.

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

The current tracked per-history comparison CSV and histogram comparison CSV
were re-read on 2026-10-20. Direct membership tests against the serialized
`emin_MeV`/`emax_MeV` intervals reproduce all 150 PHITS bin counts; no value
lands on a boundary. Maximum corrected-history absolute difference is
`7.815970093361102e-14 MeV`, absolute sum of per-history differences is
`8.223907665971808e-13 MeV`, and the total-sum absolute difference is
`1.6342482922482304e-13 MeV`. The ordinary PHITS histogram file and internal
trace are local validation inputs and are not included in Git; the tracked
CSVs are derived results. Do not infer individual histories from the ordinary
histogram alone.

The same-seed GNU baseline and diagnostic runs also had byte-identical
Dumpall and `[t-deposit]` output, supporting that the diagnostic logging hook
did not change transport for this run. PHITS's stock distribution output has
no history IDs; the internal trace supplied the IDs for the exact check.

Reproduction scripts are `validation/run_physical_reconstruction.py` and
`validation/compare_internal_history_totals.py`. The edge-based audit of the
committed derived CSVs is reproducible from the repository root with:

```powershell
python tes_cpp/dump2event/validation/audit_independent_seed_bins.py `
  tes_cpp/dump2event/validation/independent_seed_1000_history_compare.csv `
  tes_cpp/dump2event/validation/independent_seed_1000_histogram_compare.csv `
  tes_cpp/dump2event/validation/independent_seed_1000_ncol10_boundary_audit.csv `
  --source-normalization 1000 --output bin_audit.json
```

The audit was run on 2026-10-20 and returned current corrected PHITS-vs-bin
L1 `0`, corrected Dumpall-vs-bin L1 `0`, pre-correction L1 `6`, and the three
history crossings listed above. This script relies on the already-produced
PHITS-internal per-history trace export; it does not claim to regenerate the
PHITS transport or the raw diagnostic trace. PHITS diagnostic source,
binary, raw Dumpall and stopping table are local research inputs and are not
included in this branch.
