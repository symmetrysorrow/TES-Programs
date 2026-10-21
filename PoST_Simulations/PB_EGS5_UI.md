# Pb EGS5物理depositの対話UI（実験用）

`PoST_Simulations/PoST_Simulation.py` の `questionary` メニューで次を選択します。

**Pb EGS5物理depositをHDF5 v2へ変換（実験用・パルス未接続）**

画面で順番に、正常終了済みの無損失Dumpall、物理設定JSON、privateなPb HATCH停止能CSV、出力先HDF5 v2を選択します。物理設定JSONの `phits_input` が参照するPHITS入力を確認し、最後に変換を明示承認します。

通常の `tes_cpp.dump2event(..., physics="pb-egs5", physics_config=..., stopping_table=...)` を呼び出し、物理エンジンをUI側へ複製しません。終了時にはHDF5 v2として再読込し、event/particle/track/deposit件数とprovenanceの一部を表示します。全イベントのスコア完了フラグと未解決0件を確認し、巨大deposit配列全体をメモリへ読み込みません。

## 安全性

- **実験用opt-in**：デフォルトの`Dump2Event`やエネルギー補正方式は変えません。
- **対象条件限定**：PHITS 3.341、Pb region 1000 / material 2 / EGS5 medium 2、電子・陽電子cutoff 1 keV、単位重み、検証済みPb HATCH table。詳細な判定は既存converterが実施します。
- **名前を分離**：出力先初期値は`event_pb_egs5_v2.h5`。既存出力を上書きせず、旧`event.h5`と混同しません。
- **パルスへ未接続**：`Pulse_Ms()`へHDF5 v2を渡しません。両者や従来`Dump2Event`を同じメニュー選択で組み合わせると停止します。
- **物理の保証ではない**：HDF5 v2の形式・スコア完全性検査は、全吸収やSEGMENT内の熱配分が物理的に正しいことを証明しません。
- **データ保護**：配布PHITS、元Dumpall、停止能CSV、内部診断値を変更・公開しません。入力・設定不一致や未対応物理は既存CLIのfail-closedを使用します。

## テスト

`python -m pytest tes_cpp/dump2event/tests/test_pb_egs5_menu_ui.py -q`

対話処理のmock試験は、意図した物理APIへの委譲、明示確認、既存出力拒否、PHITS入力参照、HDF5形式・未解決件数を対象とします。大容量PHITS実測データをCIへ持ち込まず、科学的な全加算照合は既存のローカル検証資料を参照してください。

## mainへの採用について

研究で確認したスコア再構成コードを**実験用の選択肢として取り込む**ことと、すべてのPHITS入力への科学的な適用を認めることは異なります。mainに含めても、`pb-egs5`をデフォルトには切り替えず、旧`Pulse_Ms`がv2を誤読しない制約を維持してください。
