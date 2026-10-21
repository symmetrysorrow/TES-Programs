# 材料・停止能テーブルの仕様

実装が受け取るCSV headerは次の11列を厳密に要求する。

```csv
medium,medname,lelke,eke0,eke1,ededx0,ededx1,pdedx0,pdedx1,phits_density,egs_density
```

| 項目 | 意味 |
|---|---|
| medium / medname | PHITS material medium番号とEGS5媒体名 |
| lelke | PHITS EGS5 log-energy lookup row |
| eke0 / eke1 | index = trunc(eke1·ln(E)+eke0) のtable parameters |
| ededx0 / ededx1 | electron `dedx=ededx1·ln(E)+ededx0` |
| pdedx0 / pdedx1 | positron `dedx=pdedx1·ln(E)+pdedx0` |
| phits_density / egs_density | density provenance; validated Pb sliceでは双方9.78 |

当runの設定ではPHITS material ID 2 / EGS5 medium 2 / PB / region 1000を確認。`emin(12)=emin(13)=emin(14)=0.001 MeV`。期待medium以外のCSV行を無視し、同一mediumのmedname不一致、duplicate index、欠落lookup index、density不一致は拒否する。PDedx係数がない/インデックスが欠ける場合は値を補完しない。

ローカルHATCH係数表SHA-256 `20fd86eb0b57483bd19f20a8388d2c1ecb0431a02fa99f0933e2c3c4748514e9`。このtable自体はローカル検証由来の非公開データでありGitへ含めない。利用側はPHITS 3.341 HATCHから同じmaterial/mediumで係数を取得し、出所・SHA-256をprovenanceに保存する必要がある。新HDF5 v2にはtable hashとsource input hashがprovenanceとして書かれる。

検証対象はPB medium 2のみであり、任意材質への汎用性は主張しない。
