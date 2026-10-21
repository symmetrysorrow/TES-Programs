# 10,000履歴検証の性能と出力

- PHITS isolated diagnostic run: 2026-10-21 11:57:40–12:05:25 JST、wall time約7分45秒、exit=0。1,983,461 Dumpall records、1,954,742 tally `PRE` rows。Internal traceは約1.64 GBの非公開一時データ。
- MSVC 19.51 Release build: `tes_cpp_dump2event`, `dump2event_physical_reconstruct`, `test_physical_deposit_engine` target成功。
- Raw physical ledger scan: 2,783,824,000 byte Dumpallをstream parse、10,000 histories、1,954,742 additions、未解決0。実行時間はビルド+1 scanのログを参照。
- HDF5 v2 generation: 188.0169秒。出力は10,000 events、49,327 particles、1,973,461 tracks、1,954,742 deposits、587,741,058 bytes、unresolved=0、negative=0。
- HDF5 histogram verify: 150 bins、532 positive events、count L1=0。PHITS表示丸めとのL1 `1.554312234475219e-14`。
- Memory peakは取得していない。Raw track CSVとintermediate ledgerは`private/`に保存し、Git管理対象から除外する。
