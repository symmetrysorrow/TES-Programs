# Phase24 mortar constraint projection audit

| case | total constraint rows | TES–Stycast | Stycast–substrate | unknown |
|---|---:|---:|---:|---:|
| historical | 5560 | 1671 | 1671 | 2218 |
| Phase24 | 186 | 139 | 47 | 0 |

| interface | metric | historical mean | Phase24 mean | Phase24 / historical |
|---|---|---:|---:|---:|
| TES_to_Stycast | coefficient_count | 1.934949132e+01 | 1.790647482e+01 | 9.254235433e-01 |
| TES_to_Stycast | support_node_count | 1.934949132e+01 | 1.790647482e+01 | 9.254235433e-01 |
| TES_to_Stycast | coefficient_l2 | 8.625648212e-11 | 1.051549484e-09 | 1.219096186e+01 |
| TES_to_Stycast | support_x_span_m | 3.547513577e-05 | 1.263341113e-04 | 3.561201628e+00 |
| TES_to_Stycast | support_y_span_m | 3.530060098e-05 | 1.173368125e-04 | 3.323932432e+00 |
| Stycast_to_substrate | coefficient_count | 1.942190305e+01 | 1.593617021e+01 | 8.205256802e-01 |
| Stycast_to_substrate | support_node_count | 1.942190305e+01 | 1.593617021e+01 | 8.205256802e-01 |
| Stycast_to_substrate | coefficient_l2 | 8.602593684e-11 | 3.158699681e-09 | 3.671799224e+01 |
| Stycast_to_substrate | support_x_span_m | 3.543116404e-05 | 2.061671469e-04 | 5.818808171e+00 |
| Stycast_to_substrate | support_y_span_m | 3.524666792e-05 | 2.409943725e-04 | 6.837366103e+00 |

Constraint rows are classified from their primal support nodes. The raw constraint-row count is not itself a conductance; it diagnoses mortar tessellation and projection stencil changes. Phase24 has fewer but much larger support patches; the total coefficient L1/L2 mass is reported in the JSON to distinguish area aggregation from an actual projection-strength change.
