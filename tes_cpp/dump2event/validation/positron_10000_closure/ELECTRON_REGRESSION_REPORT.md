# 電子系回帰検証

電子計算式は変更せず、PHITS 3.341 Pb EGS5 table loaderの共通構造にPDedx列を保持する項目を追加した。既存electron coefficient lookup、NCOL=10 boundary、NCOL=11 cutoff、dead electron secondary、strict exact-zero branchを維持した。

同一runのPHITS internal additionと今回のraw-only ledger:

| 分岐 | 加算数 | 内部最大絶対差 (MeV) |
|---|---:|---:|
| EGS5 electron continuous | 1,901,831 | 1.1449174941446927e-16 |
| EGS5 boundary | 83 | 9.8879238130678004e-17 |
| EGS5 cutoff | 38,266 | 5.865533753146579e-15 |
| dead electron secondary | 1,310 | 0 |
| strict exact-zero | 11,448 | 0 |

電子の既存100履歴回帰: 38,077/38,077加算、最大差5.7976626971489864e-15 MeV、加算ID・history・particle・region差0、checked positive position差0。MSVC C++ physical-engine unit test成功。`pytest tes_cpp/tests`の結果は`PUSH_AND_CI_REPORT.md`に記録する。
