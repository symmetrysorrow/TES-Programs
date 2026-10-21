# 明示的に未検証の物理・設定

物理モードはresearch-onlyかつ明示的な`--physics pb-egs5 --hdf5-version 2` opt-inである。以下は、現時点の1,000履歴/2 MeV/Pb/単位重みの内部加算joinによって一般化されない。

- 非単位重み、履歴途中のweight change、splitting/roulette/importanceの寄与。
- Pb以外、region 1000以外、Dumpall material 2以外、EGS5 medium 2以外、cutoffが1 keV以外。
- PHITS 3.341以外の版、`[t-deposit] output=dose`、空間meshなど異なるscore/position semantics。
- arbitrary `esmin`, `eletb(1)`。閾値はPHITS 3.341の入力/既定と実行ソースから解決され、未知ならfail-closedする。
- 2 MeV datasetに現れなかった、または内部加算と再構成を個別分類していない輸送・score branches。
- 全入射エネルギー収支、脱出・残留エネルギー、物質内全吸収状態。HDF5 `complete`はスコア対象再構成の完了であり、全吸収確認ではない。
- 履歴重みが変わる場合の`output=deposit`履歴分布を唯一の正解にすること。
- PR #6の既定`E_deposit`式への適用。本研究モードの結果を未検証範囲の代替式やゼロ埋めに使用しない。

現時点で、2 MeV/1000 runの既存CLI HDF5 v2はその同じrunの全5,978,912内部加算と一致した。ただし独立seed 10,000 runは未完了のため、それに固有の処理分類はまだ結果なし。
