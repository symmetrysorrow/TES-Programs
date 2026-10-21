# Pb EGS5物理deposit UI（dump2event唯一の変換方式）

`PoST_Simulations/PoST_Simulation.py --interactive` の
**「Pb EGS5物理depositをHDF5 v2へ変換（パルス未接続）」** を選択します。

旧親粒子エネルギー差方式はmainのCLI・Python API・メニューから削除済みです。
物理方式選択メニューはなく、Dump2Eventは必ず検証済みPb/EGS5物理再構成を呼び出します。

順番に無損失Dumpall、物理設定JSON（`phits_input`でPHITS入力を参照）、
private Pb HATCH停止能CSV、出力先HDF5 v2を選択します。
初期出力名は`event_pb_egs5_v2.h5`です。既存ファイルは上書きしません。

変換前の明示確認、入力ファイルの存在、条件外・未対応物理のfail-closed、
正常終了後のHDF5形式・スコア完了・未解決件数確認を実施します。
内部物理式はUIに複製せず、同じ`tes_cpp.dump2event` Python APIを呼び出します。

**HDF5 v2→Pulse_Msはまだ未接続です。** SEGMENTに蓄積したエネルギーをTESへ
配分するモデルは今後`event2pulse`で別途検証してください。
同時実行による誤接続も拒否します。全吸収や輸送全体の完全性を意味しません。

テスト：`python -m pytest tes_cpp/dump2event/tests/test_pb_egs5_menu_ui.py -q`
