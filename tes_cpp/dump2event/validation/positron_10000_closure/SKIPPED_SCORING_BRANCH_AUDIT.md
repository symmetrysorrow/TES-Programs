# 未実装・今回の検証slice外の分岐

完全なinternal tally joinでは、今回のPb/単位重み/同一10,000履歴で全1,954,742加算を再構成器と照合した。これは値と計上位置のslice全体比較であり、分岐ごとの物理式すべてを別seed/設定で一般化した検証ではない。

明示的に未検証または未対応の範囲:

- splitting, roulette, importance等で履歴中に重みが変化する場合。特にPHITS event-distributionの重み保持規則を一般化していない。
- QCLUSTS(7)が保守的guard `1e-6 MeV` 以下のdead charged secondary。
- Pb region 1000/material 2/EGS5 medium 2以外。
- 光子反応別・原子緩和別・領域脱出別の個別のDumpall-only物理ledger。全タリーadditionsの内側照合は通ったが、それら各過程のraw-only分類までは完了していない。
- 動的物質・ユーザー定義 `dedxfnc`、異なる tally unit、未取得/欠落したHATCH stopping table。
- positive continuous-deposit intervalの中でPHITSが選んだ代表点。raw endpoint intervalは保存するが、中点代表位置とは仮定しない。
- 独立seedでのpositronを含む全過程の internal-addition一致。

未対応値を0にしたり、clip/normalize/再配分したりしない。現在の検証は既知の明示設定を前提に、対象runのscore additions全件を特定条件で照合するもの。
