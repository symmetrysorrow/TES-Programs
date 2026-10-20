# Rest-mass field unit audit

Existing raw-only analysis of the separate official Windows 2 MeV, 1,000-history Dumpall reports the dead positron archival-child `QCLUSTS(5)` value as `0.510998902 GeV` per child. This is `510.998902 MeV`, close to the electron rest mass; it is not a 0.510998902 MeV value. Do not add it as local deposit: it is a child rest-mass record and may coexist with transported annihilation photons. The generated live photon children each carry their own energy in their relevant QCLUSTS energy fields.

Because the GNU samples do not parse, this unit audit was not repeated for the current 20-history runs. Source units and historical raw values are documented as provenance, not a new run mapping. No rest-mass contribution was added to the experimental converter.
