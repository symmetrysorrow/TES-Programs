# NOT_SCORED、exact-zero、UNRESOLVEDの意味

| 状態 | PHITS内部加算 | HDF5 deposit | 今回の判定基準 |
|---|---:|---:|---|
| `NOT_SCORED_BY_PHITS` | なし | 作らない | PHITS 3.341のscore条件を満たさず関数に到達しないとソースと実行設定から証明 |
| exact-zero scored addition | あり、tlv=0 | 0 MeV行を保持 | 内部TALLYSEQ/PRE-POST実記録に存在 |
| `UNRESOLVED` | 不明または再構成不能 | 完成HDF5として公開しない | スコア対象か値・位置を証明できない |
| `INVALID` | 該当せず | 失敗 | 不正構造、非有限値、設定不一致 |
| positive scored addition | あり、tlv>0 | 再構成値を保持 | ソース対応分岐と同一runのaddition IDを照合 |

NOT_SCORED候補はaudit CSV/provenance件数で追跡する。deposit数に加えない。実在するゼロ加算とは別に数える。2 MeV/1000では全5,978,912 PHITS加算の一対一対応を確認し、NOT_SCORED候補1、actual zero additions 4,797、unresolved 0、invalid 0だった。候補はNOCAS 398の子index 0で、QCLUSTS(7)=7.5933446597176157e-7 MeV < eletb(1)=9.9999999999999995e-7 MeV。親反応に関連する別加算は親レコードから生じたもので、この子候補の加算ではない。10,000履歴の値は同run join終了後に追記する。
