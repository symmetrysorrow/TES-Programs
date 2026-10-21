# Gamma-process coverage

The current engine scores charged-particle transport and supported dead charged secondaries using the verified Pb EGS5 branches. Photon interactions are not converted directly into local deposits. Their energy is expected to be transported by charged secondaries or to escape, but each relevant PHITS score branch must be checked against same-run internal additions.

| Condition | Status in this integration task |
|---|---|
| Co-60 1.173/1.332 MeV photons, Pb region 1000 | Historical tally exists; new same-run raw/internal proof pending provenance repair |
| Photoelectric / Compton / Rayleigh | Not independently revalidated through the integrated CLI |
| Pair creation and e+/e- transport | Historical e+ branches exist; independent seed pending |
| Positron annihilation photons and subsequent transport | No 1.022 MeV local rest-mass deposit is added by this mode; gamma score coverage not complete |
| Photon escape and Pb/vacuum boundary | Not independently validated through this CLI |
| 511 keV and 2 MeV / thick and thin Pb fixtures | Existing separate diagnostic assets are not treated as integrated CLI proof |

Any unsupported score candidate causes the complete-output path to fail rather than silently use the old parent-minus-children formula.
