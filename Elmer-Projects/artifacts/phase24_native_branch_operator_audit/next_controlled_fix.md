# Next controlled diagnostic fix (proposal only; not executed)

## Target

Refine the **TES/Membrane conformal trace at z = 192 µm** (`TES__zmin` ↔ `Membrane_SiNx__zmax`, BID 1104/1305) in the TES-footprint edge band, s = max(|x|, |y − 1000 µm|) ∈ [230, 270] µm, to at most 12 µm. The historical trace is 11.7 µm with 174 nodes on the edge line. Best Phase24 is currently 38.6–39.0 µm with 48 nodes.

Keep everything below z = 191 µm unchanged:

* Membrane_Si1
* SiO2_1
* frame
* Si_2 / SiO2_2
* bath faces

Those planes are already 10.7–12.7 µm in the band, which is historical-equivalent.

## Why this one

* **Operator split.** 81% of the +4.337% gap comes from the isothermal TES-edge loading term. The same network loaded by uniform footprint flux differs by only 0.82%.
* **Mesh.** z = 192 is the only plane in the TES-edge band that is still 3.3× coarser than historical.
* **Zone dissipation.** The resistance missing in best is concentrated in s = 230–270 µm, about 170% of the net ΔR, in the sheets fed through this trace.
* **Membrane_SiNx sheet.** Its +15% stiffness, about 0.9% of G, comes from joining a 39 µm top to a 10.7 µm bottom across 1 µm. A finer top face should reduce this bulk term too.

## Expected result (prediction, not measured)

* Removing the edge-loading term (1.0348 → ~1.00) gives G/G_hist ≈ 1.008.
* Leaving 20% of that term in place gives G/G_hist ≈ 1.015.

Both are inside the 2–3% band.

## Gate for accepting the variant

* **Must stay unchanged:** node/face counts of every body and plane below z = 191 µm, the bath (SiO2_2 zmin) faces, and the frame.
* **Allowed to change:**
  * the TES mesh (its native G_dissipation is 1.6e-12 W/K);
  * the TES-top/Stycast mortar and the Stycast/abs mesh. They form a dead end with net Q ≈ 1e-18 W and G_dissipation ≈ 3e-17, so they cannot move steady G by more than about 1e-9 relative.
  * the Membrane_SiNx interior.
* **Must be re-audited:** after the three frozen-power CPU/MUMPS points, run `analyze_phase24_native_branch_operator_audit.py`. The edge-loading factor should approach historical, and sheet A's NtD error should fall below about 3%.

## Only if |G_variant/G_hist − 1| ≤ 0.03

Run the full nonlinear CPU/MUMPS steady validation. HYPRE/GPU stays NO-GO until thermal parity is confirmed.
