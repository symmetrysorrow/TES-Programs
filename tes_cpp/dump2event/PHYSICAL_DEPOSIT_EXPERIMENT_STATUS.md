# Experimental PHITS deposit reconstruction status

This branch contains a separate opt-in C++ ledger executable. It does not
change the default `dump2event` deposit formula. It is a research experiment,
not a complete PHITS deposit converter and is not ready for PR #6.

## Validated slice

- PHITS 3.341, Pb region 1000 / material and EGS5 medium 2.
- 1 keV electron transport cutoff and unit transport weights only.
- PHITS-derived EGS5 stopping table, supplied as a local input rather than
  committed data.
- EGS5 electron continuous loss (NCOL 13/14, JCOLL 14), electron boundary
  transport (NCOL 10, ITYP 12), electron cutoff (NCOL 11, ITYP 12), dead
  electron release above a conservative scoring guard, and exact `E == EC`
  zero scores.

The NCOL 10 branch was added after tracing `tdepstregEVENT` and `egs5edxde`:
the PHITS boundary tally scores `ein-ecc`, and the EGS5 path gives
`ecc=ein-dedx(E)*step`. The eight positive boundary terms omitted by the first
raw-only version explain the prior independent-seed histogram discrepancy.
Details and the row-level audit are in `INDEPENDENT_SEED_ROOT_CAUSE.md` and
`validation/independent_seed_1000_ncol10_boundary_audit.csv`.

The previous L1=6 consisted of three adjacent-bin moves in NOCAS 246, 365, and
461, caused by those omitted NCOL=10 terms. Exact PHITS serialized histogram
edges and the Fortran half-open comparisons rule out a boundary-rounding
explanation for those cases. After correction, the tracked PHITS-internal
comparison has maximum per-history absolute error `7.815970093361102e-14 MeV`
and histogram L1=0. The separately quoted `2.22e-15 MeV` is the raw-only
implementation-to-implementation maximum, not the PHITS-internal maximum.

## Independent-seed result

On the existing 1,000-history independent-seed case, the revised raw-Dumpall
ledger was compared to PHITS 3.341's per-addition internal diagnostic trace.
There were 201,971 additions on each side, no unresolved rows, and no
per-history difference above `1e-12 MeV`. Total absolute difference was
`1.6342482922482304e-13 MeV`; maximum per-history difference was
`7.815970093361102e-14 MeV`. Both histories and reconstructed deposits used
the same `[0,1.5)` range and 0.01 MeV half-open bins; histogram L1 difference
was zero. See `validation/independent_seed_1000_history_compare.json`.

This is one PHITS 3.341 configuration and seed. It validates these additions
and the selected material/weight/cutoff case only. It does not validate
arbitrary materials, non-unit weights, positron annihilation, photon or other
scoring branches, spatial voxel splitting, or all transport boundaries.

## Regression baseline retained

The existing 100-history proof artifact reports 38,077/38,077 tally additions
matched, maximum error `5.7976626971489864e-15 MeV`, zero addition/history/
region mismatches, and zero positive-deposit coordinate mismatches. The
10,000-history structural parser read 1,983,461 records through NCOL=2. The
current MSVC 19.51 build succeeds, the C++ physical-engine test passes, and
the Python suite passes 23/23. The default converter regression processed
1,000 histories successfully. These checks do not expand the physics scope.

## Failure-closed behavior and limitations

The program prints `status=EXPERIMENTAL_PARTIAL` even if every row in a sample
is covered. Missing physics settings are not inferred. A dead secondary with
energy `<= 1e-6 MeV` is unresolved because PHITS' exact `esmin`, which affects
the tally threshold `min(esmin,1e-6 MeV)`, is not in Dumpall. Unsupported
NCOL=11 particles are emitted as unresolved even if their raw `E` and `EC`
are equal. Other unsupported branches and weights are not normalized or
clipped. The prototype writes a ledger, not event JSON/HDF5.

The 100-history 38,077-addition result and the independent-seed comparison are
local validation artifacts; the large PHITS data, PHITS source, diagnostic
builds, and nonpublic stopping table are excluded from Git.
