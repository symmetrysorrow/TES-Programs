# Positron and annihilation accounting

The 0.511 MeV mono-gamma run produced zero positrons, so it provides no pair-production or annihilation validation. In the 2 MeV, 100-history run, HDF5 readback contained 13 distinct `ITYP=13` positron tracks and the raw export contained 3,624 positron transport rows. The physical ledger classified 3,624 positron-continuous additions totaling 3.267614070158497 MeV and 13 dead-positron-secondary additions totaling 0.01130702759570578 MeV; all 3,637 entries had positions.

The ledger keeps positron kinetic transport/local scoring separate from rest mass and from annihilation photons. A positron's 511 keV rest-mass term plus the target electron rest mass is represented through the annihilation products and subsequent transport/absorption, not as an additional 1.022 MeV local deposit at the vertex. This avoids double-counting emitted photons. The 100-history internal addition join found no energy, position, or region mismatch for the ledger's scored additions.

**Limit:** 13 tracks are a small sample. This does not prove every flight-in-annihilation, in-flight annihilation, at-rest annihilation, photon escape, or downstream photon interaction branch. The 1,000-history results are pending in the run report and any absent branch will remain unverified.
