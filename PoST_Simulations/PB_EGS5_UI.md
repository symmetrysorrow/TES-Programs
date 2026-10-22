# Pb EGS5 → HDF5 v2 → Pulse_Ms

PoSTの対話メニューでは`Dump2Event`がPHITS 3.341/Pb/EGS5の物理
deposit再構成によりHDF5 v2を出力します。従来の親粒子エネルギー差変換
および新しいv1出力はありません。検証済みのPb条件以外は変換時に拒否します。

変換にはDumpall、同じrunに対応するPHITS入力・実行ログ、検証済みPb HATCH
停止能テーブルが必要です。ソースファイルは通常同じディレクトリから
自動関連付けされ、bundledの停止能表はSHA-256を確認して使用します。
明示的なprivate CSVの指定も可能です。変換先の既存ファイルは上書きしません。

## 物理イベントをパルスへ変換する

`Pulse_Ms`は旧HDF5 v1（点付与）と新HDF5 v2（物理付与）を判別します。
v2では`POINT`を`score_xyz_cm`に、`SEGMENT`を
`start_xyz_cm`〜`end_xyz_cm`の線分に割り当てます。線分が複数の
absorber blockに跨る場合、各blockでの通過長に比例してエネルギーを
分配します。総エネルギーは保存します。x座標の左右反転モードも
始点・終点の双方に適用します。

**これは線分内で一様dE/dsを仮定する近似モデルです。** PHITSの
スコア再構成検証だけでは線分内の詳細な熱密度分布は確定できません。
現状では各depositの異なる時刻も解決できないため、全blockの参照波形を
同じ時間起点で線形重畳します。正しい熱輸送の検証と区別してください。

不完全/未解決のイベント、Pb以外のmaterial、領域1000以外の付与、
単位/座標系不整合、absorber外の位置については、黙ってエネルギーを捨てず
失敗します。既存のpulse出力はエラー時には保持されます。
`all_absorption=undetermined`はパルス生成の妨げにはなりませんが、
全吸収であるという意味ではありません。

入力イベントの配置は以下に対応します。

- 各位置ディレクトリ: `<event-root>/<posi>/event.h5`
  または `<event-root>/<posi>/event_pb_egs5_v2.h5`
- 位置が1つだけなら `--event-root=/path/to/any-v2-file.h5` を直接指定可能
- 対話メニューでDump2EventとPulse_Msを同時選択する場合は、
  `input.json`の`position`が1位置のとき変換結果を直接受け渡します
- 複数位置では位置ごとに物理変換した後にPulse_Msを実行してください

`pulses.h5`は既存の`posi2pulse`による単一block参照波形です。
今回の変更ではその計算モデルには手を加えず、eventの空間配分を変更しました。
結果は各`<output>/<posi>/pulse_MS.h5`に保存され、
`source_event_format_version`、`energy_allocation_model`、
`deposit_time_model`などの属性を含みます。

テスト: `python -m pytest PoST_Simulations/tests/test_event2pulse_v2.py -q`

## Bi research status

Biは現段階で物理再構成が検証された通常出力としては未承認です。
PbのHATCH係数やモデルをBiへ代用しないでください。
詳細は`tes_cpp/dump2event/BI_VALIDATION_PLAN.md`を参照してください。
