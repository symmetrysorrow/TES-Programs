# Full absorption status

All histories in all three HDF5 files carry `all_absorption=undetermined`. A complete supported deposit subset does not prove full energy balance or absence of particle escape. The schema's structural `complete` field is true for 100/100 and 1,000/1,000 because those files contain no unresolved rows for the configured scoring slice; for 10,000 histories it is true for 9,991 and false for 9 histories containing unresolved entries. This does not mean full transport or all-absorption was confirmed.

No event is marked `confirmed`. Leakage and energy-balance evidence across unsupported particle branches is insufficient for a confirmed/rejected decision.
