# HDF5 v2完了・収支フィールドの意味

HDF5 v2の`/events/complete`は、対象とするスコア候補を処理し、`UNRESOLVED`を残さずHDF5を閉じられることを表す。物理的な全吸収や入射エネルギーの全収支が証明されたことは意味しない。

2 MeV / 1000履歴のreadback実測:

- `/events/complete`: true 1000/1000
- `/events/unresolved_count`: 0が1000/1000
- `/events/all_absorption`: `undetermined` 1000/1000
- provenance `transport_physics_completeness`: `undetermined`
- provenance `energy_balance_status`: `not assessed`
- provenance `complete_semantics`: no in-scope scored candidate remains `UNRESOLVED`; proven not-scored candidates are excluded from deposits and preserved in audit

そのため、このrunは出力・スコア対象加算の対応が完了している一方、dump2eventで輸送全過程の物理エネルギー収支や全吸収状態を確定したものではない。10,000履歴についても同じ基準でflag分布を調べ、全吸収を根拠なく`confirmed`へ変更しない。
