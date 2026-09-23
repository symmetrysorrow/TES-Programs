# Phase24 mortar Schur/control audit

## D-block check

- The captured saddle systems have an exactly zero lower-right D block.
- Therefore the literal B^T D^-1 B expression is undefined.
- The valid constraint Schur complement after eliminating primal x is -B K^-1 B^T.

| case | center constraints | T(0.95P) mK | T(1.00P) mK | T(1.05P) mK | G derivative W/K | G / historical |
|---|---:|---:|---:|---:|---:|---:|
| historical | 5560 | 166.200114399 | 167.052750880 | 167.905387362 | 1.878294461631e-08 | 1.000000000 |
| original_phase24 | 186 | 160.967076184 | 161.544289413 | 162.121502641 | 2.774542062577e-08 | 1.477160328 |
| phase24_interface_refined | 1721 | 164.998259546 | 165.787638324 | 166.577017101 | 2.028813576885e-08 | 1.080136059 |

## Controlled interface refinement

- Original Phase24 G: 2.774542062577e-08 W/K
- TES/Stycast interface-refined G: 2.028813576885e-08 W/K
- Reduction from original Phase24: 0.268775340
- Refined/historical G ratio: 1.080136059

The refinement changes only the local TES top / first 1 um Stycast mesh in this existing diagnostic case. The outer Stycast-substrate interface remains coarse. The large G reduction therefore strongly implicates the TES-Stycast interface discretization or its local FE/mortar coupling, although it is not a pure algebra-only mortar test because the local element mesh is refined as well.
