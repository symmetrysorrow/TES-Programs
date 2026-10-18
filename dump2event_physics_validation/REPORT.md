# Dumpall → dump2event 負エネルギー付与の物理検証

## 結論

**最も強い個別の説明は、陽電子消滅で生成される約1.022 MeVを、親陽電子の運動エネルギーだけから差し引いていることです。確度はこの反応例について高いものの、負値12,626件の主因全体を説明する確度は低いです。** 最大負値1.020998 MeVの行は、親陽電子の運動エネルギー1 keVに対して、二つの光子の運動エネルギーが各0.510998902 MeVで、第三の停止陽電子が0.001 MeVです。これは消滅で解放される静止質量エネルギー約1.022 MeVと一致します。

残る負値の大半は電子（ITYP=12）のEGS5反応（NCOL=14, JCOLL=14）で、子粒子はほぼ全て電子2個です。CSVの計算式では、記録された親エネルギーより二次電子の運動エネルギー合計が大きくなっています。ただし、ここから直ちにPHITSの輸送が物理的にエネルギーを作ったとは言えません。親子の粒子重み、E/ECの項目対応、束縛電子・原子反跳のエネルギー、EGS5の記録規約を生Dumpallと同一反応で照合する必要があります。

## 実行データと対象範囲

- 入力: `H:\hata\1332_142_136\1332keV_6\dumpall.dat`、ASCII Dumpall、116,820,667,876 bytes。
- PHITS出力 `phits.out`: PHITS 3.341、`negs=1`、`e-mode=2`、`dumpall=-1`。`trace.inp` は1.173 MeVと1.332 MeVのCo-60光子を各強度1で指定しています。
- 監査CSVは10,000履歴の処理結果で、全反応監査1,909,651行、負値12,626行。元CSVは `D:\temp\analysis\pr6-reaction-audit-10000\` にあり、このフォルダにはコピーしていません。
- PRソース参照: commit `36dec9f139918d233faf0876736f43d1bbdb965c` の `tes_cpp/dump2event/src/converter.cpp`。この実装の補正前局所付与値は `deposited_before_collision - transported_secondary_energy_sum` です。
- Dumpall公式マニュアルはPHITS 3.34版を参照しました。NCOL=13/14でNCLSTSと二次粒子情報が出力され、JCLUSTS(4)が状態、QCLUSTS(7)が運動エネルギー(MeV)、QCLUSTS(8)が重みです。ここでは**監査CSVに重み列がなく**、重みを考慮した再集計は未実施です。

## 集計結果

| 指標 | 結果 |
|---|---:|
| 負値反応行 | 12,626 |
| 負値反応を1件以上含む履歴 | 533 |
| 絶対値の中央値 / 90パーセンタイル | 1.739 keV / 22.552 keV |
| 最大絶対値 | 1.020998 MeV |
| `parent_E - parent_EC` の最大絶対値 | 1.0×10⁻¹⁵ MeV |
| 親ITYP=12, NCOL=14, JCOLL=14 | 12,589件（負値全体の99.7%） |
| 負値反応で最初の二次粒子だけでも親エネルギー超過 | 5,608件 |
| 負値行中、二次粒子にdead状態あり / 閾値未満あり | 5件 / 0件 |

全反応監査CSVを分母にした負値率は、ITYP=12, NCOL=14, JCOLL=14で0.660%、ITYP=13, NCOL=14, JCOLL=14で2.373%、ITYP=14, NCOL=14, JCOLL=13で0.574%でした。親エネルギー帯別の値は [negative_deposit_statistics.csv](negative_deposit_statistics.csv) にあります。

負値の二次粒子構成は、親電子→電子+電子が12,456件、親電子→電子+光子が133件、親陽電子→電子+陽電子が29件、親光子→電子+光子が3件です。陽電子→陽電子+光子+光子の5行はdead陽電子を含み、後述の大きな例とは別の構成です。

完全保持された75履歴は、1.173 MeVが34件、1.332 MeVが41件で、各グループ全件が `energy_consistent=1` でした。履歴合計の最大残差はそれぞれ1.71×10⁻¹³ MeV、1.88×10⁻¹³ MeVです。したがって、先に疑われた「34履歴の1.173 MeV取り違え」は履歴レベルでは説明になりません。個々の負値反応が存在しても、全履歴の収支は閉じています。

## 代表例と物理的解釈

代表行とCSVに記録された値は [reaction_examples.csv](reaction_examples.csv) にまとめました。

1. **履歴318209・粒子2・反応7（電子EGS5）**: 親1.006438 MeV、二次電子1.047554 MeVと0.001381 MeV、合計1.048935 MeV、局所付与−0.042497 MeV。最初の二次粒子単独で親エネルギーを超えます。これは記録エネルギー・重みの扱いを調べる優先例です。
2. **履歴119953・粒子3・反応404（陽電子消滅候補）**: 親陽電子0.001 MeV、光子0.510998902 MeVが2個、dead陽電子0.001 MeV。実装はdead陽電子のエネルギーを差し引かず、輸送光子1.021998 MeVを親の1 keVから差し引くため、−1.020998 MeVになります。光子2個の合計は電子・陽電子の静止質量エネルギー2mₑc²と一致します。PHITS公式マニュアルではJCOLL=2はparticle decay、QCLUSTS(7)は運動エネルギーと定義されています。
3. **親電子→電子+光子**: 負値133件。光子分岐、電子の連続損失、重みの扱いを個別に見る必要があります。
4. **親陽電子→電子+陽電子**: 負値29件。粒子種と運動エネルギーだけでは、対生成・消滅に関係する静止質量項を閉じられない可能性があります。
5. **非負値比較**: 典型的な電子EGS5反応を `reaction_examples.csv` に含めました。平均差だけでなく、同じ親エネルギー、同じ反応種、重み条件で比較する必要があります。

### 生レコード照合の制限

PHITSの公式定義では、NCOL=4行にNOCAS・NOBCH・RCASC・RSOUINがあり、NOCASとRCASCは別の値になり得ます。履歴119953は生行5464050–5464051でNOCAS=RCASC=119953を確認できました。同一履歴内の生行5475209–5475240にJCOLL=2、NCLSTS=3、0.510998902 MeV光子2個とdead陽電子0.001 MeVがあり、監査CSVの消滅候補と娘粒子・エネルギーが一致します。ただし粒子番号3と履歴内反応番号404は、生行から独立に数え直していません。

電子反応の代表例、event318209/particle2/reaction7は、監査CSV上で二次電子1.047554 MeVと0.001381 MeVです。該当エネルギー値のraw行を一意に見つけられず、line-numbered raw matchは未完了です。したがってこの電子例は監査CSVの報告値であり、生Dumpallから独立確認済みとは扱っていません。どちらの例も**履歴・反応の完全な生レコード照合には未完了項目が残ります**。

## `[t-deposit]` 照合とPHITS再実行

既存の `trace.inp` には `[t-deposit]` がなく、計算フォルダにも履歴別deposit tally出力を確認できませんでした。そのため `tally_comparison.csv` は作成していません。PHITS実行ファイルも環境内で確認できず、複製入力での再実行は未実施です。PHITSマニュアル上、`[t-deposit]` の `output=deposit` とdump出力は履歴別の領域付与エネルギーを記録できますが、元計算の各反応点と同じ粒子ごとの二次記録をそのまま出す機能ではありません。再実行で比較するなら、セル別・履歴別の総depositと、履歴の完全収支をPR監査の集計に照合する設計が必要です。

## 問いへの回答

1. **二次粒子の運動エネルギー合計が親の記録エネルギーを超えるか**: はい。監査CSVでは12,626件すべてがそう計算され、5,608件は最初の二次粒子単独でも親を超えています。ただし、重みや親子間の物理的エネルギー項を含めた完全な収支は未検証です。
2. **停止粒子の扱いミスか**: 負値全体の説明にはなりません。dead二次粒子がある負値行は5件だけで、残り12,621件は監査上輸送扱いです。大きな陽電子例ではdead陽電子0.001 MeVを局所付与へ加えない扱いも寄与しますが、支配項は二つの511 keV光子です。
3. **正負の付与が打ち消し合い、75履歴の収支が閉じるか**: はい。75件すべてが源線別の履歴エネルギー判定を通過し、1.173 MeV群34/34、1.332 MeV群41/41です。ただし、この閉合はPR実装の同じDumpall解釈による合計で、独立なPHITS tally検証ではありません。
4. **追加調査の優先順位**: ①監査CSVのイベント番号をNOCAS/RCASCと対応付けて対象反応の生レコードを再抽出、②QCLUSTS(8)と親WTを含む重み付き二次エネルギー合計、③E/ECおよびJCLUSTS状態の列対応、④陽電子消滅の静止質量エネルギー、⑤同条件の履歴別`t-deposit`総量、の順です。
5. **修正方針**: [PROPOSED_FIX.md](PROPOSED_FIX.md)。この時点では既存 `E_deposit` の仕様変更を提案しません。まず独立tallyと重み込みの収支を確立し、陽電子消滅で媒質電子の静止質量をどの収支項に置くか決めます。

## 参照資料

- [PHITS User’s Manual Ver. 3.34 (JAEA PDF), §5.2.25 Dumpall](https://phits.jaea.go.jp/manual/manualE-phits334.pdf) — NCOL、NOCAS/RCASC、E/EC、NCLSTS、JCLUSTS、QCLUSTSの定義。
- [PHITS公式マニュアル §6.7.22 Dump function](https://phits.jaea.go.jp/manual/PHITS-en/chapters/tally-format.html) — `[t-deposit] output=deposit` の履歴別dump説明。
- [PHITS公式マニュアル §5.2.13 EGS5](https://phits.jaea.go.jp/manual/PHITS-en/chapters/sections/parameters/parameters_13.html) — EGS5設定と電子・陽電子・光子のcutoff設定。


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

?????????????? [reaction_index_mapping.csv](reaction_index_mapping.csv) ??????
