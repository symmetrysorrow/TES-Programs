# 最終判定

## 結果

今回新たに実施した同一seed診断で、従来の342 unresolved行は、326 EGS5 positron NCOL14輸送、7 annihilation直前のNCOL13 positron輸送、9 dead positron secondaryに分類できた。別に旧HDF5がゼロ扱いした1,460件の`E==EC` positron transport additionsを特定した。全1,802件の加算はDumpallと明示的Pb HATCH PDedx係数から再現し、同一run internal tlvと1e-12 MeV tolerance内で一致した。

修正実験版のHDF5 v2は同じraw 10k Dumpallから新規作成。1,954,742/1,954,742加算がinternal traceと順序照合され、unresolved/negative/value/region/positive-position mismatchは0。Per-history histogramはPHITS 150-bin tallyとcount単位で一致し、L1=0。PHITSの印字丸めに対する連続値の残差L1は1.55e-14。

## 実装可否

新しい明示的な実験用物理再構成器とHDF5 v2 builderを、PHITS 3.341 Pb region 1000、EGS5 medium 2、1 keV cutoff、単位重みというscope内の研究用実装としてpushする根拠が得られた。標準`dump2event` CLIとPR #6の反応差分式は変更しない。

これは物理全般の完成ではない。動的weights、他材料・領域、energy cutoff以下dead charged particles、光子反応ごとの独立raw-only台帳、独立seedでの全加算内部照合、連続lossの内部代表点は未検証。PHITS全吸収とHDF5 `complete`も同じ意味とはしない。PR #6への統合は別の独立検証が必要。

## 安全策

- 値のclip、history normalization、エネルギー再配分なし。
- 重み・粒子種・region/material・PHITS version条件を限定し、外れ値はUNRESOLVEDに保持。
- 標準converterのdefault動作は変更しない。物理器は別の明示的実験用 executable/pipeline。
- PHITS source, coefficient table, 10k Dumpall, diagnostic trace, generated HDF5 and raw-level CSV are excluded from Git.
