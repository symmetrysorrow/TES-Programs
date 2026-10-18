from pathlib import Path
p=Path('dump2event_physics_validation/REPORT.md')
section=r'''

## 2026-10-19 追記: 電子EGS5反応の生Dumpall照合

今回の再照合では、監査CSVの `event_id` はNCOL=4ヘッダーの第1列 `NOCAS`、`reaction_index` は履歴全体で数えたNCOL=13/14反応の通番として扱った。PR #6の `converter.cpp` はNCOL=4の `column[0]` を `case_number` に代入し、その値を `event_id` に出力している。RCASCをevent_idとして検索する方法は誤りだった。

### 電子EGS5: NOCAS=318209、粒子2、第7反応

実DumpallのNOCAS=318209ヘッダーは生行612–613にあり、値は `NOCAS=318209, NOBCH=1, RCASC=318209, RSOUIN=12`。同一履歴のNCOL=13/14レコードを順に数えると第7反応は生行849のNCOL=14、親粒子ID=2、ITYP=12、JCOLL=14である。生行849–884がこの反応と2つの子粒子の完全なレコードを含む。

| 生レコード項目 | 生Dumpall値 | 監査CSV値 | 差 |
|---|---:|---:|---:|
| 親E (MeV) | 1.00643799692798 | 1.00643799692798 | 0 |
| 親EC (MeV) | 1.00643799692798 | 1.00643799692798 | 0 |
| 親WT / OLDWT | 1 / 1 | 監査CSVに列なし | — |
| 子0: ITYP / JCLUSTS(4) / QCLUSTS(7) / QCLUSTS(8) | 12 / 0 / 1.0475541999824101 / 1 | 12 / 0 / 1.0475541999824101 / 重み列なし | エネルギー差0 |
| 子1: ITYP / JCLUSTS(4) / QCLUSTS(7) / QCLUSTS(8) | 12 / 0 / 0.0013812126696086671 / 1 | 12 / 0 / 0.0013812126696086671 / 重み列なし | エネルギー差0 |

親EまたはECと子の運動エネルギー合計の非加重差は `1.00643799692798 − (1.0475541999824101 + 0.0013812126696086671) = −0.0424974157240388 MeV`。重み付き差も、`親E×OLDWT − Σ(QCLUSTS(7)×QCLUSTS(8))` および `親EC×WT − 同じ子合計` のどちらも−0.0424974157240388 MeVである。各重みは1なのでこの反応では両収支は一致する。監査CSVの局所付与−0.042497415724038756 MeVとも丸め誤差内で一致する。

### 陽電子消滅: NOCAS=119953、粒子3、第404反応

この履歴のNCOL=4は生行5464052–5464053で、`NOCAS=119953, NOBCH=1, RCASC=119953, RSOUIN=1102`。全NCOL=13/14レコードの通番404は生行5475195のNCOL=13、親粒子ID=3、ITYP=13、JCOLL=2。生行5475195–5475240で親レコードと3つの二次レコードを確認した。

親のE=0.001000000000001 MeV、EC=0.001 MeV、WT=1、OLDWT=1。二次粒子は光子2個（各 `ITYP=14, JCLUSTS(4)=0, QCLUSTS(7)=0.510998902 MeV, QCLUSTS(8)=1`）とdead陽電子1個（`ITYP=13, JCLUSTS(4)=-1, QCLUSTS(7)=0.001 MeV, QCLUSTS(8)=1`）。したがって全子粒子を含む非加重エネルギー合計と重み付き合計はいずれも1.022997804 MeVで、親Eとの差は−1.021997803999 MeV。dead状態でもQCLUSTS(8)は1である。

監査CSVの従来の「輸送二次粒子だけ」の合計は、dead陽電子を除いた光子2個の1.021997804 MeVである。これを親Eと比較すると−1.020997803999 MeVになり、監査CSVの−1.020997803999999 MeVと丸め誤差内で一致する。つまり、その局所付与値は全二次粒子を含む収支ではなく、dead二次粒子を除いた輸送粒子の収支に対応する。

### NOCASとRCASCが異なるヘッダー、およびPHITS出力条件

先にRCASC=318209で拾った生行3と64はそれぞれ `(NOCAS=191679, NOBCH=1, RCASC=318209, RSOUIN=2)`、`(NOCAS=156535, NOBCH=1, RCASC=318209, RSOUIN=1)` で、監査イベント318209のNOCASではない。正しいイベントヘッダーは生行613である。今回対象ヘッダーまでの生行1–5,515,858の範囲で、RCASC=318209かつNOBCH=1の別NOCASヘッダーを計11件確認した。全行は [observed_nocas_rcasc_mismatches.csv](observed_nocas_rcasc_mismatches.csv) と [history_id_mapping.csv](history_id_mapping.csv) に記録した。この件数は対象ヘッダーまでの走査範囲であり、116.8 GB全体の不一致件数を表すものではない。

元の `phits.out` はPHITS 3.341、`maxcas=1,000,000`、`maxbch=1`、`negs=1`、`e-mode=2`、`dumpall=-1`、`icntl=0` を記録する。元の `trace.inp` は光子源の `e-type=8` に1.173 MeVと1.332 MeVを各強度1で指定している。PHITS 3.34マニュアルのNCOL=4定義では `RCASC = NOCAS + maxcas×(NOBCH−1)`。同式は現行3.37公式マニュアルにも掲載されている。maxcas=1,000,000、NOBCH=1ならRCASCはNOCASと一致するはずなので、生行3・64などの値とは整合しない。

確認できた出力仕様、`maxcas`、`maxbch`、入力条件から、この不一致を説明する条件は特定できなかった。PHITS 3.341固有の別定義、ファイル連結・再起動・並列出力の影響のいずれかを示す証拠も今回の資料では得られていない。**NOBCH=1でNOCAS≠RCASCとなる理由は未解決**として記録する。RCASCの一致をイベントID対応の根拠に使わない。

### 追加成果物

- [raw_dumpall_pr6_electron_annihilation_excerpts.txt](raw_dumpall_pr6_electron_annihilation_excerpts.txt): 元Dumpallの行番号付きヘッダーおよび対象反応ブロック。
- [history_id_mapping.csv](history_id_mapping.csv): 監査event_id、NCOL=4のNOCAS、NOBCH、RCASC、行番号の対応。
- [parent_child_energy_weight_comparison.csv](parent_child_energy_weight_comparison.csv): 親E/EC/WT/OLDWT、各子粒子の状態・運動エネルギー・重み、加重/非加重収支と監査CSVの比較。
- `validation_scripts/scan_target_headers.py`, `inspect_reaction_headers.py`, `write_pr6_raw_artifacts.py`: 抽出と反応通番の確認手順。

この追記ではdeposit計算式を変更していない。電子例は生レコード完全一致を確認したが、その結果が示すのは記録値と現在の差分計算の一致であり、物質・束縛電子を含む最終的な物理エネルギー会計の説明までは確定しない。
'''
with p.open('a',encoding='utf-8',newline='\n') as f:f.write(section)
