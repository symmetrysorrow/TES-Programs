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

## Bi research status

The material is read from the PHITS physics contract JSON (`region_material`).
Pb is the previously validated single-material route. Bi is recognized but
rejected **before** any HDF5 conversion because material-specific HATCH
coefficients and same-run PHITS score-addition validation have not been
established. The UI must never silently substitute Pb settings for Bi.
See `tes_cpp/dump2event/BI_VALIDATION_PLAN.md`.

## Pb停止能表は原則として自動選択

現在の対話UIはPb停止能CSVを毎回選択しません。
`tes_cpp/dump2event/data/pb_egs5_phits3341.csv`に正確な実測済み
HATCH係数を入れて配布・インストールすれば、自動でSHA-256を検証して使用します。
現時点では元CSVの公開再配布条件が未確認で、公開Gitには本体がありません。
よって必要な場合はCLIの`--stopping-table`、Python APIの
`stopping_table=...`で私有表を指定してください。
