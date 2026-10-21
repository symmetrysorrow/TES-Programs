# NOCAS 398 dead-secondaryの原因

状態: **スコア対象外の判定は確定**。同一親DumpallレコードにPHITS内部加算があることと、当該dead secondary index=0に独立加算が作られないことを、2 MeV/1000の全件joinと同run TALLYSEQで区別して確認済み。

## raw Dumpallから得た値

通常CLIの読み取り対象は同一runの `D:\temp\gamma_pb_cli_physics_20261021\seq_n1000_2mev\gamma_2mev_thickPb_dumpall.dat`。同CLIが出力した監査候補CSVの唯一の行:

```text
NOCAS=398, raw physical line=39594769, parent particle ID=180, parent ITYP=14,
NCOL=13, JCOLL=13, secondary index=0, ITYP=12, JCLUSTS(4)=-1,
JCLUSTS(5)=-1, QCLUSTS(7)=7.5933446597176157e-7 MeV,
QCLUSTS(8)=1, region=1000, material=2,
xyz=(-0.89993733682112143,-0.0080077241081570032,-3.4727350725342681) cm
```

監査候補CSVが記録する当該分類は `NOT_SCORED_BY_PHITS`。変換器はこの候補にaddition IDもHDF5 deposit行も与えない。

## 判定根拠

実測 `eletb(1)=9.9999999999999995e-7 MeV` に対して
`QCLUSTS(7)=7.5933446597176157e-7 MeV < eletb(1)`。
PHITS 3.341 `talls08.f:4929-4930` はdeadかつcharged二次粒子について、さらに厳密な `QCLUSTS(7)>eletb(1)` を満たしてからscorer呼出へ進む。この粒子はその条件を満たさず、score関数へ到達しない。

しきい値超過量は `Q-eletb1 = -2.4066553402823838e-7 MeV`。単に内部加算が見つからなかったため除外したのではなく、ソース分岐と同run runtime値で不到達を証明した。

## 内部trace照合

同一runの `internal_additions.csv` をstream検索すると、同じ親raw source lineの内部PRE/POST/TALLYSEQ additionは2件だった。両方ともNOCAS 398, particle ID 180, NCOL=13, JCOLL=13で、増分は `0` と `0.0010959999999999859 MeV`。2つとも親輸送区間の同じ座標端点を持つ。これは同じ親反応行に複数のadditionがあることを示すが、対象dead child index=0 (0.759 keV) の独立depositではない。score gateがこの子を除外するソース/runtime証拠と、実加算が2件ある事実は両立する。子index=0をHDF5加算IDに誤って結びつけない。

通常CLIが出力したHDF5と同一run内部加算5,978,912件をjoinした結果、加算IDの欠落・重複0、座標/region/particle/event不一致0だった。従って親レコードにdepositが存在することだけを、特定dead childのQCLUSTS(7)がdepositされた証拠として扱わない。

内部加算全体とHDF5の最大値差は4.680023729664029e-14 MeV。タリー条件の実測はstrict `QCLUSTS(7)>eletb(1)`。候補差分は `-2.4066553402823838e-7 MeV`。対象子の加算IDは存在しないのが正しいPHITS実行分岐の結果であり、未解決0件と同一の意味ではない。
