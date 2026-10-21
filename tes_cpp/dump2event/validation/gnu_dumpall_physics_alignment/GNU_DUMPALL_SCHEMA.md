# PHITS 3.341 Dumpall structural schema

Source of truth: installed PHITS 3.341 `C:\PHITS\src\ovly15.f`, routine `dumpdat` (logical write order) and `cpcasdump` (thread-file merge). ASCII Dumpall is list-directed numeric text; physical line lengths are not a safe field schema.

| Logical group | Values | Present when |
|---|---:|---|
| NCOL marker | 1 | file start |
| History header | NOCAS, NOBCH, RCASC, RSOUIN | NCOL=4 |
| Particle prefix | particle ID, material ID, ITYP, KTYP, JTYP, MTYP, RTYP, OLDWT | NCOL 4–15 |
| QS | 1 | ITYP 12 or 13 |
| Region/lattice levels | region before/after, lattice levels before/after | every data record |
| Lattice arrays | 5 values per level, before then after | level > 0 |
| Geometry | COSTHA, angular direction 1–3, NSURF | every data record |
| Particle counters | name and 3 counters | every data record |
| Transport | WT and direction cosines U,V,W | every data record |
| State at E point | E,T,X,Y,Z | every data record |
| State at EC point | EC,TC,XC,YC,ZC | every data record |
| Spin and state | SPX,SPY,SPZ,NZST | every data record |
| Collision block | NCLSTS; if positive: MATHZ,MATHN,JCOLL,KCOLL and each child’s ICLUSTS, JCLUSTS(0:8), QCLUSTS(0:12), JCOUNT(1:3) | NCOL 13/14 |
| Terminal markers | NCOL=3 and NCOL=2 | run/batch and end markers as emitted |

For the particle prefix, the GNU sample must contain eight complete numeric values. `0.000` at the end of an 80-character record is syntactically numeric but cannot establish whether it is a truncated RTYP or a complete zero-valued field. The following region record cannot be shifted into this tuple without corrupting every subsequent field. The parser therefore validates group counts and fails closed.
