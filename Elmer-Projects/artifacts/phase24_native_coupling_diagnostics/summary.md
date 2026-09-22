# Phase24 native TES–membrane coupling diagnostics

- Isolated mesh: `mesh_phase24_native_coupling_diag`; the 179 TES–membrane contact nodes were duplicated on the membrane side, so the native mortar interface is nonconforming and has zero shared nodes.
- Physics, materials, TES law, circuit constants, bath temperature, and CPU MUMPS backend were retained. Only the diagnostic TES power was frozen at 0.95, 1.00, and 1.05 times the historical checkpoint power.
- Native TES–membrane rows: 179 in every case; native full saddle constraints: 365.
- Symmetric bracket: ΔT = 1.154420524 mK, ΔQ = 3.203005460e-11 W, G_eff = 2.774556926e-08 W/K.
- The CSV contains the actual native Cλ reaction integrals and temperatures for all three solves.
