# Refined Phase24 thermal-resistance decomposition

This read-only analysis reuses the existing historical and refined h=10 um mortar native MUMPS captures. No solver, production mesh, material, or TES law was changed.

The endpoint budget is TES -> Membrane_SiNx -> bath-side SiO2_2 -> bath. The middle term includes the intervening SiO2/Si/SiNx/Membrane stack and any parallel network paths; it is not a unique single-interface resistance.

- TES_to_Membrane_SiNx: historical `1.562892228e+07 K/W`, refined mortar `1.380961473e+07 K/W`, delta `-1.819307555e+06 K/W`, share `0.458`
- Membrane_SiNx_to_bath_side_SiO2_2_including_stack: historical `3.761077765e+07 K/W`, refined mortar `3.545379679e+07 K/W`, delta `-2.156980862e+06 K/W`, share `0.542`
- bath_side_SiO2_2_to_bath: historical `1.156211688e+02 K/W`, refined mortar `1.156299668e+02 K/W`, delta `8.798062449e-03 K/W`, share `-0.000`
- total_TES_to_bath: historical `5.323981556e+07 K/W`, refined mortar `4.926352715e+07 K/W`, delta `-3.976288408e+06 K/W`, share `1.000`

Per-body temperature slopes are in `body_temperature_slopes.csv`; the three-point temperature profiles are in `body_temperature_profile.csv`.
The TES->Membrane and bath-side substrate->bath endpoint terms are separated directly. Any dominant residual in the middle term points to the membrane/substrate network realization, not TES-Stycast mortar formulation.
