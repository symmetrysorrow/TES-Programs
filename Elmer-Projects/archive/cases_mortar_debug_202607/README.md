# Mortar-debug case archive (2026-07)

モルタル境界の不安定性を切り分けるために使った実験ケース群。結論は
`docs/README_TES_elmer.md` に反映済みで、現行ワークフローでは使わない。

- `case_constant_power_fixed*.sif` — merged メッシュ上の熱伝導率固定・部位除去実験
  (メッシュ: `mesh_shifted_merged`)
- `case_constant_power_fixed_refined.sif` — 細分化テストメッシュ
  (メッシュ: `legacy/meshes/mesh_refined_test` へ移動済み)
- `case_constant_power_unmerged*.sif` — unmerged メッシュ + モルタル段階的有効化
  (メッシュ: `legacy/meshes/mesh_shifted_unmerged` へ移動済み)

再実行する場合はリポジトリルートにコピーし、`Mesh DB` のパスを現在の
メッシュ位置に合わせて修正すること。
