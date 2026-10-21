# Pb/Bi material support — implementation and validation gates

## Status (PHITS 3.341; 2026-10-22)

| Material | Strict contract | EGS5 ledger | Real same-run PHITS addition comparison | Finished HDF5 v2 |
|---|---|---|---|---|
| Pb | Implemented (existing region 1000, material 2, medium 2, cell density 9.78) | Validated subset | **Yes**, 2-MeV/1,000-history 5,978,912 additions | Allowed with strict contract |
| Bi | Implemented and testable with *synthetic* coefficient fixtures | **Research only** | **Not performed** | **Rejected** until independently validated |

The existing C++ stopping-power equations (EDedx/PDedx) are shared; this does
**not** mean the Bi material-dependent coefficients or PHITS score branches are
validated. The synthetic Bi unit-test numbers have no physical significance.

## Configuration

The existing Pb physics contract is unchanged, including its canonical resolved
configuration hash. Pb may optionally specify \`cell_density_g_cm3: 9.78\`.

For Bi, the JSON needs the standard PHITS 3.341/1-keV/unit-weight/scoring
fields plus the following **explicit** values. **No real numeric values are
inferred or provided**:

- \`region_material: "Bi"\`
- \`region_id\`: integer obtained from the actual PHITS input
- \`dumpall_material_id\`: corresponding positive material ID
- \`egs5_medium\`: the actual HATCH medium index, **not necessarily** equal to the material ID
- \`cell_density_g_cm3\`: positive value taken from the actual PHITS source input
- \`stopping_table_sha256\`: exact SHA-256 of the matching private HATCH CSV
- \`phits_input\`: path to the corresponding single-material PHITS input
- \`esmin_MeV\` and \`dead_secondary_score_threshold_MeV\`: sourced from PHITS 3.341 runtime settings

The source validator compares \`mat[material_id] Bi 1\`, the configured scoring
cell's density and material, and the \`[t-deposit]\` scoring region, while
retaining the existing 1-keV cutoffs, unit weights, and 150-bin tally definition.

Missing density, inconsistent material/region/medium, other elements/alloys,
unrecognized weights, or nonmatching PHITS input are rejected. A Bi contract
that is otherwise syntactically sound is still rejected by the **normal HDF5
driver**, with a clear "not yet PHITS-oracle validated" message.

## Research-only ledger (NO published Bi HDF5)

The C++ research helper can read a private Bi HATCH table and independent
material-ID/medium mappings. Its explicit interface is:

\`\`\`text
dump2event_physical_reconstruct <complete_bi.dumpall> <private_bi_hatch.csv>
    <candidate_ledger.csv> <non_scored_audit.csv> <resolved_eletb1_MeV>
    Bi <region_id> <dumpall_material_id> <egs5_medium> --research-bi
\`\`\`

Its output is **unvalidated candidate data**, not an endorsed physical deposit
dataset. It does not use the PHITS internal trace as an input. All formulas,
threshold branches, energy indices, and score positions need independent
verification against PHITS internal \`[t-deposit]\` additions.

## Next validation without an expensive 10,000-history run

1. Create an **isolated** PHITS 3.341 Bi input and matching EGS5 HATCH table;
   log density, composition, region, material ID, medium, cutoffs and executable
   hashes. Do not modify \`C:\PHITS\`.
2. Start with **100 histories** and obtain lossless Dumpall plus the internal
   score addition trace from the **same run** using a nonperturbing diagnostic
   build. Check complete NCOL=2 terminator and header arity.
3. Independently compare every occurrence (history, particle, region, process,
   addition order), including electrons, positrons, zero score and
   \`NOT_SCORED_BY_PHITS\` threshold candidates. Record maximum per-addition
   errors, coordinates, POINT and SEGMENT endpoints.
4. Verify history totals and 150-bin \`[t-deposit]\` bins *in addition* to
   one-to-one internal score matching. Never fit the distribution or fabricate
   energy. Keep all absorption \`undetermined\` where not proven.
5. Resolve all candidate branches, then repeat at a modest independent seed
   before proposing a **separate, reviewable change** to permit Bi HDF5.
   No change to the Pb validation contract, no large run required.
6. Handle mixed Pb+Bi transport and steps crossing material boundaries only
   after a separate multi-region score/segment validation. Current HDF5 v2
   contains a region ID per deposit and a single-material provenance contract;
   it is **not** a verified mixed-material heat allocation format.

## Compatibility

\`dump2event\` continues to produce HDF5 v2 only. The old signed-difference
converter is not restored. Pb works with the existing private Pb HATCH table;
Bi cannot fall back to Pb data. The PoST interactive menu identifies the
material in the JSON and rejects Bi **before** conversion until oracle closure
is demonstrated. \`Pulse_Ms\` remains disconnected from v2; event2pulse will
be handled separately.
