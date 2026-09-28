# Phase24 mortar-neighborhood element stiffness audit

- Stycast conductivity used: 2.69094e-06 W/(m K)
- k*A/h is a local parent-element normal stiffness scale, not the assembled mortar conductance.
- Temperature gradients come from saved 1.00P native full_x_after captures mapped through the result Temperature Perm table.

| interface side | faces historical | faces Phase24 | area ratio | mean h ratio | sum kA/h ratio | sum exact face-block K ratio | mean gradient ratio | heatflow proxy ratio |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| TES_to_Stycast:side_a | 4274 | 4048 | 1.000000000 | 0.444049483 | 2.440949557 | 1.179321118 | 15.955317563 | 17.006602007 |
| TES_to_Stycast:side_b | 3214 | 3177 | 1.000031039 | 0.423447074 | 2.613625453 | 1.102256573 | 0.548934397 | 0.541747072 |
| Stycast_to_substrate:side_a | 3214 | 74 | 0.980221734 | 0.381395469 | 2.824412588 | 0.206289025 | 0.529794831 | 0.538085200 |
| Stycast_to_substrate:side_b | 11444 | 314 | 1.000000000 | 5.146337829 | 0.205071681 | 0.139877563 | 0.178142168 | 0.290544626 |

The face-area and local-stiffness comparisons are split by side because the mortar interfaces are nonconforming. A large side-to-side gradient change with near-equal area would point to local element geometry or temperature interpolation; a large area or kA/h change would point to face/height integration. The exact face-block K metric is computed from the element conduction matrix and is the stronger comparison.
