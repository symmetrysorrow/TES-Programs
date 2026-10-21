# Negative, zero, duplicate, and missing-score audit

- The 10,000-history candidate HDF5 reported 0 negative and 0 unresolved rows; previous same-seed internal comparisons reported 0 mismatches above `1e-12 MeV` and 0 region mismatches. The provenance issue in this integration task prevents those figures from being promoted to independent same-run confirmation.
- A regression guards the prior positron `E == EC` false-zero case by ensuring positron transport routes through PDedx before the exact-zero branch.
- Unit-weight rejection is enforced for all parsed parent and secondary weights. The complete CLI rejects unresolved candidates; it does not clip, normalize or fallback.
- Addition count classes sum to 1,954,742 in the stored candidate ledger. Two zero-chord positron rows are reclassified from the old exact-zero count; their ledger zero is known, while their direct internal-tlv association remains unverified.
- The existing C++ HDF5 writer checks schema/event append integrity; there is not yet a full integrated raw-to-HDF5 duplicated addition-ID audit on a new independent run.
