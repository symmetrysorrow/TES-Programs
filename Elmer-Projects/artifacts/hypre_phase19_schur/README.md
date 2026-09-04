# HYPRE Phase19 Schur diagnostic artifacts

These JSON files are the checkout-visible snapshots used by the Phase19
documentation.  They are independent one-rank CPU diagnostics on the saved
explicit mortar matrix:

- `exact_schur.json`: SuperLU K factorization and dense Schur correctness oracle.
- `diag_schur.json`: `diag(K)^-1` Schur approximation.
- `ilu_proxy.json`: sparse ILU K-action proxy.
- `mgr_best.json`: the targeted MGR best result copied from the handoff status.

The JSON values are not production acceptance claims.  The exact dense Schur
is an oracle, and the MGR result reached its iteration limit above the gate.
The matrix-free Elmer action validator is
`scripts/analysis/validate_matrix_free_schur.py`; it requires the four
`*_vN.dat`/`*_yN.dat` files emitted by an Elmer diagnostic run.
The current harness attempt is recorded in `matrix_free_diagnostic_status.json`;
it emitted only the first pair before the surrounding one-step solve failed.
