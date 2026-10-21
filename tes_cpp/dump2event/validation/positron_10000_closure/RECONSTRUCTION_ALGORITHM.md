# Pb / EGS5 陽電子のraw-only再構成

## 通常陽電子輸送寄与

対象はPHITS 3.341 EGS5のITP=13 primary trackで、NCOL=13または14。PHITSが選ぶHATCH energy rowは

```text
ek = max(E_raw, emin(13))
logE = ln(ek)
index = trunc(eke1 * logE + eke0)
PDedx = pdedx1[index] * logE + pdedx0[index]
L = sqrt((xc-x)^2 + (yc-y)^2 + (zc-z)^2)
tlv = OLDWT * PDedx * L
```

となる。`emin(13)`、medium/material、係数は明示入力で与える。未解決行をこの式へ無条件に流さず、region=1000/Pb、medium=2、PHITS 3.341、unit weights、およびNCOL/ITYPを確認してから使う。内部式はPHITS 3.341の`src/talls08.f`と`src/egs5edxde.f`で確認した。これは生Dumpallと既知のPHITS物理設定からの独立再計算であり、内部trace値を計算入力には使用しない。

## Dead陽電子secondary

対象はNCOL 13/14のdead charged secondaryで、`ITYP=13`, `JCLUSTS(4)<0`, charge `JCLUSTS(5)!=0`, finite `QCLUSTS(7)>1e-6 MeV`, unit QCLUSTS(8)の検証済みsliceのみ。

```text
E_local = QCLUSTS(7) - 0
position = (QCLUSTS(10), QCLUSTS(11), QCLUSTS(12))
```

これは停止後の残余**運動**エネルギーであり、静止質量を含まない。ソースはQ7を`ein`、0を`ecc`に置き、normal branchで`tlv=ein-ecc`とする。実装は保守的なDumpall-independent guard以下の値を再構成せずUNRESOLVEDに残す。

## 位置と領域

連続損失は衝突点の一点depositにせず、Dumpallの始点・終点で示される輸送区間として保存する。内部診断の区間端点とraw端点を比較する。dead secondaryはraw child座標を点として保存する。別途PHITSスコア代表点と同一と確認できない区間の中点は作らない。region 1000から外れる値をregion 1000へ押し込まない。

## 重みと残りの会計

この実装のスコープはunit weightのPb検証sliceである。連続輸送式はPHITS内部のOLDWTを使用し、dead-secondary通常枝はQCLUSTS(8)を掛けない。動的重みの履歴分布は本実装の保証外。消滅で生じる511 keV光子は新規transport粒子であり、dead陽電子のQ7から推測して同位置depositに加えない。
