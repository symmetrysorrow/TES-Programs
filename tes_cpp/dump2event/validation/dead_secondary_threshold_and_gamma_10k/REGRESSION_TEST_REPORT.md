# 回帰試験記録

現ブランチの物理設定・閾値変更後、2026-10-21にMSVC 19.51 Releaseビルドと回帰テストを再実行した。

## ビルド

```powershell
cmake --build build\physdeposit-vs --config Release --target dump2event dump2event_physical_reconstruct test_physical_deposit_engine -j 4
```

結果: **成功**。3 targetがReleaseでビルドされた。

## 試験

```powershell
ctest --test-dir build\physdeposit-vs -C Release --output-on-failure
python -m pytest tes_cpp\dump2event\tests -q
```

- CTest: 1/1 passed。
- pytest: 34 passed, 10 skipped。
- スキップ10件は専用の実データ/PHITS診断アーティファクト依存。ローカル一時runを作成しないと実行されないため、passとして数えていない。

再ビルド後にも同じコマンドを実行し、同じ結果を得た。既定CLI実データ回帰は、以前実行済みの511 keV/100, 511 keV/1000, 2 MeV/100および今回の2 MeV/1000のPHITS/HDF5照合を併記する。通常物理モードは明示オプションでのみ呼ばれ、コードの通常JSONモードをこの変更では変更していない。

なお、PHITS 2 MeV/10,000履歴計算は現在稼働中で、その後の通常CLI、全加算joinは未実施。これらを回帰成功数に含めない。
