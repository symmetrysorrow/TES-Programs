# Gamma, electron and positron process classification

## 0.511 MeV, 100 histories

This source is below the pair-production threshold. The raw transport export contains 250,805 records, all mapped to unique Dumpall record sequence IDs. HDF5 contains 4,955 distinct particles: 4,954 electrons and one source photon; positrons: 0. The physical ledger has 250,938 additions: 250,590 positive and 348 exact zero; negative, unresolved, and unsupported additions: 0. This run validates photon interaction/deposit handling and secondary-electron transport, not positron behavior.

## 2 MeV, 100 histories

Distinct transported particles in HDF5: 355 photons, 11,025 electrons, 13 positrons (11,393 total). Raw transport rows by `ITYP`: photons 657, electrons 545,293, positrons 3,624. The difference between raw rows and distinct particles reflects the multiple transport records per track.

| Parent ITYP | Ledger process | Additions | Sum (MeV, unweighted) | Position available |
|---:|---|---:|---:|---:|
| 12 | EGS5 electron continuous | 534,262 | 83.14071710840383 | 534,262 |
| 12 | electron cutoff | 11,011 | 11.076523887160953 | 11,011 |
| 12 | electron boundary | 10 | 0.02143825466993763 | 10 |
| 12 | dead electron secondary | 4 | 0.0012814053952309168 | 4 |
| 13 | EGS5 positron continuous | 3,624 | 3.267614070158497 | 3,624 |
| 13 | dead positron secondary | 13 | 0.01130702759570578 | 13 |
| 14 | dead electron secondary | 270 | 0.7899957718337252 | 270 |
| 14 | exact zero | 445 | 0 | 0 |

The parent-`ITYP=14` dead-electron-secondary rows arise from electron products scored at photon records; they are classified by the scoring parent and are not photon energy-difference estimates. The ledger records positron kinetic transport and dead secondary contributions separately; annihilation photons remain transport energy until a later local interaction or escape. Neither the 1.022 MeV rest-mass term nor the annihilation photons are inserted as an extra local deposit at the annihilation vertex.

These are observed counts from one 100-history 2 MeV sample. They do not establish every positron termination or annihilation branch.

## 2 MeV, 1,000 histories

The raw transport export contains 5,977,725 track records: photons `ITYP=14`: 6,841; electrons `ITYP=12`: 5,937,907; positrons `ITYP=13`: 32,977. The raw-only physical ledger contains 5,978,912 reconstructed additions and one unresolved `unsupported_secondary_threshold` child. The ordinary CLI rejected HDF5 publication because of that candidate. The per-process addition and energy totals are listed below after streaming aggregation; they describe candidate local scoring from the same raw run, not an HDF5 output.

| Parent ITYP | Ledger process | Count | Reconstructed sum (MeV) |
|---:|---|---:|---:|
| 12 | EGS5 electron continuous | 5,818,560 | 897.5017379593426 |
| 12 | EGS5 electron cutoff | 119,229 | 119.9401250448232 |
| 12 | EGS5 electron boundary | 59 | 0.2132197669313079 |
| 12 | dead electron secondary | 93 | 0.33744762803479433 |
| 13 | EGS5 positron continuous | 32,977 | 31.078712889268864 |
| 13 | dead positron secondary | 118 | 0.1085182561453185 |
| 14 | dead electron secondary | 3,101 | 10.414865777894413 |
| 14 | exact zero | 4,775 | 0 |
| 14 | unresolved secondary below threshold guard | 1 | blank, excluded from sum |

The raw-only supported ledger sum is 1,059.594627322437 MeV; the PHITS internal `tlv` sum is 1,059.5946273224363 MeV. These sums agree within `6.82e-11 MeV`, but the single unresolved candidate remains visible in the ledger and blocks HDF5 publication. Rest-mass energy is not added at annihilation; the positron and annihilation-photon energy rows remain separated by the observed later transport/interactions.
