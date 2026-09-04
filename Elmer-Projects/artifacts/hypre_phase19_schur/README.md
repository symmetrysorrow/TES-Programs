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
`matrix_free_validation_current.json` is the current four-vector result and
is an explicit failure snapshot: vector parity passes, but the strict action
gate does not (`max action relative error = 5.0716313e-7`).
`superlu_parity_cpu.json` additionally records K backward residuals,
componentwise backward errors, `B^T v`/`B K^-1 B^T v`/`Dv` stage parity, and
raw SHA-256 fingerprints for K/B/Bt/D.  The diagnostic emits the corresponding
`*_K|B|Bt|D.triplets` and `*_vN|btN|kuN|bkuN|dvN|yN.dat` files in the project
root.  The raw block fingerprints mismatch for K/B/Bt; only D matches.
Lower/full CPU, GPU, MPI, and transient promotion tests were intentionally not
run after this failed gate.
