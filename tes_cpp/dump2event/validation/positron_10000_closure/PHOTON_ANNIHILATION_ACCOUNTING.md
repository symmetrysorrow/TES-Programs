# 陽電子静止質量と消滅光子の会計

- positron primaryのEGS5 stopping lossは運動エネルギーに対するPDedx式である。
- dead positron secondary scoreはQCLUSTS(7)（運動エネルギー）−0。QCLUSTS(5)で表される静止質量は含まれない。
- annihilationで発生する2本の511 keV光子はITYP=14の輸送二次粒子として残り、消滅位置のlocal depositとして加算しない。
- 光子のtransportと後続する光電・コンプトン・対生成・原子緩和分岐の個別物理台帳は今回新規に監査していない。全10kの[t-deposit]加算値は内部scoreとの照合に含まれるが、これをもって各光子過程のraw-only粒子エネルギー会計や位置モデルを一般検証済みとは呼ばない。
- 従って511 keV photonsの物質内吸収は、PHITSが後続荷電粒子輸送として実際に計上した場合に限り加算される。光子脱出はdepositへ足さない。

Rest massをdead e+ Q7へ加える処理は根拠がなく、本研究実装では行わない。
