# Physics configuration contract

The JSON contract is explicit and exact-match validated:

- `phits_version`: `3.341`
- `phits_input`: path to the PHITS input, resolved relative to the JSON file
- `region_id`: `1000`; `region_material`: `Pb`
- `dumpall_material_id`: `2`; `egs5_medium`: `2`
- electron and positron cutoffs: `0.001 MeV`
- `weight_mode`: `unit`
- `energy_unit`: `MeV`
- `deposit_mode`: `t-deposit-output-deposit`

The referenced input must explicitly define a 64-bit `bitrseed`, `icntl=0`, `maxbch=1`, `negs=1`, `e-mode=2`, pure Pb as material 2, region 1000 with material 2/density -9.78, validated electron/positron cutoffs, `totfact=1`, the 1.173 and 1.332 MeV unit-intensity source lines, and the matching Pb `[t-deposit]` region output. The seed and source-input hash are recorded in HDF5 provenance. The reconstruction executable independently rejects non-unit parent WT/OLDWT and secondary QCLUSTS(8) anywhere in the transport stream. Stopping-power medium/density/index validation remains in `Egs5StoppingTable`.

Dynamic splitting/roulette/importance weights, other materials/regions, other cutoffs, alternate tally units, and other PHITS versions are rejected or unsupported. This checks declared input/settings and observed Dumpall weights; it does not independently query a live PHITS executable's runtime version.
