# Phase24 mortar constraint projection audit

| case | total constraint rows | TES–Stycast | Stycast–substrate | unknown |
|---|---:|---:|---:|---:|
| historical | 5560 | 1671 | 1671 | 2218 |
| Phase24 | 5811 | 5766 | 45 | 0 |

| interface | metric | historical mean | Phase24 mean | Phase24 / historical |
|---|---|---:|---:|---:|
| TES_to_Stycast | coefficient_count | 1.934949132e+01 | 9.517169615e+00 | 4.918563210e-01 |
| TES_to_Stycast | support_node_count | 1.934949132e+01 | 9.517169615e+00 | 4.918563210e-01 |
| TES_to_Stycast | coefficient_l2 | 8.625648212e-11 | 1.370771103e-11 | 1.589180395e-01 |
| TES_to_Stycast | support_x_span_m | 3.547513577e-05 | 5.265781404e-05 | 1.484358351e+00 |
| TES_to_Stycast | support_y_span_m | 3.530060098e-05 | 4.830696062e-05 | 1.368445842e+00 |
| Stycast_to_substrate | coefficient_count | 1.942190305e+01 | 1.608888889e+01 | 8.283888992e-01 |
| Stycast_to_substrate | support_node_count | 1.942190305e+01 | 1.608888889e+01 | 8.283888992e-01 |
| Stycast_to_substrate | coefficient_l2 | 8.602593684e-11 | 3.287795559e-09 | 3.821865451e+01 |
| Stycast_to_substrate | support_x_span_m | 3.543116404e-05 | 2.185693042e-04 | 6.168843448e+00 |
| Stycast_to_substrate | support_y_span_m | 3.524666792e-05 | 2.368804050e-04 | 6.720646773e+00 |

Constraint rows are classified from their primal support nodes. The raw constraint-row count is not itself a conductance; it diagnoses mortar tessellation and projection stencil changes. Row count and support-patch size can move in opposite directions; the total coefficient L1/L2 mass is reported in the JSON to distinguish area aggregation from an actual projection-strength change.
