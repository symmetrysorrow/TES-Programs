# Phase24 interface primal residual accounting

- Dimensions: 96955 = 96769 primal + 186 constraints.
- Transient: BDF1 dt=1.80000000e-05 s; pulse contribution=0.
- Constraint check max |B_t x_s-g_t|: 5.85019726e-14.
- Worst primal row: 4547 (interface_primal), eta=4.00210774e-03, r=-6.64092651e-12.
- Worst interface primal row: 4547.
- Raw formula accounting error: 2.01999848e-09; closed error after adding steady restart residual: 4.13590306e-25.
- Dominant unbalanced term: steady_restart_primal_residual (L2=2.36279611e-09); raw DeltaKx and Delta b are large but cancel.
- DeltaKx L2=7.07001687e-08; DeltaMortar L2=1.64770010e-14; Delta b L2=7.07001729e-08; r L2=2.36279609e-09.
- Mass/history max mismatch: 0.00000000e+00; physical K delta*x L2=9.51951764e-15.
- Top-20 categories: {'interface_primal': 20}; top-100 categories: {'interface_primal': 99, 'bulk_primal': 1}.

Raw captures do not contain separately tagged Joule, boundary, or pulse RHS vectors. rhs_contribution_diff.csv reports the exact delta_b-history remainder; pulse is zero by input. Assembly code locations are in interface_primal_residual.json.
