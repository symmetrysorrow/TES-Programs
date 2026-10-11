# Phase24ソース再現性の調査結果

## 結論

**検証済み隔離Phase24ソースは、公式Elmer `release-26.1` commit `5a8de867068be0568f09af40fb90fee300dfbede` と、今回作成した単一の順序付き差分 `01_phase24_source_from_release_26_1.patch` から再現できました。** 新しいRelease buildも完了し、全2,237工程が成功、HeatSolveを含むsolver/moduleを生成しました。

これにより「どの上流ソースへどのパッチを当てれば、今ある検証済みPhase24 source snapshotになるか」は再現可能になりました。ただし旧ローカルrevision `9916c3f05` のGit objectは現在のリポジトリ・公式remoteから取得できず、旧累積パッチのHeatSolve基底とも一致しません。したがって歴史的なcommit ancestry自体を復元したとは報告しません。

キャッシュ既定値はOFFのままです。今回のRelease buildはソース再現の確認用で、SIFを使ったruntime数値検証や75 ms計算は実行していません。本番DLL/SIFとVector Assembly設定は変更していません。

## 履歴の対応関係

```text
ElmerCSC/elmerfem release-26.1
  commit 5a8de867068be0568f09af40fb90fee300dfbede
  ├─ 公式HeatSolve blob 1cae7fd66cd4...
  └─ 01_phase24_source_from_release_26_1.patch
       └─ 16 tracked source filesを変更
          └─ 検証済みの隔離Phase24ソース状態
             HeatSolve raw SHA-256 bde260c7...fea6193d
             HeatSolve normalized Git blob e9b6351b1ddf...
             └─ Release -O3 build成功

旧のruntime表示 Rev: 9916c3f05
  └─ 完全なGit commit object / checkoutは現在のremote・workspaceにない
     └─ 旧 elmer_native_source.patch のHeatSolve preimage bc3335c28...
        は公式5a8...のHeatSolve 1cae7fd6...と異なる
```

旧 `elmer_native_source.patch` は履歴版 `0f740c55` から復元し、SHA-256 `30153dd02451392e272dfac6f6ddbca7adb54e3c52d38ba7d4a9b51ca471afbc` を記録しました。公式5a8…への通常適用はHeatSolveの約line 666で失敗しました。HeatSolveを除外すると残りのhunkは適用できますが、その出力は現在の隔離ソースと15ファイル中12ファイルで異なります。この旧パッチを再現手順には使いません。

旧パッチのindexはHeatSolve preimage `bc3335c28` → output `b050ca5c7` を要求します。公式5a8…のpreimageは `1cae7fd66` です。対象の `9916c3f05` はテストfixtureと文書に短縮revisionとして現れますが、現在のElmer cloneには存在せず、公式GitHub remoteにもrefとして広告されず、`git fetch ... 9916c3f05` も解決しませんでした。

リポジトリに追跡された `source/HeatSolve.F90` はGit blob `603c847f0de44b2a88d8801655b599d26c9ddf22`（SHA-256 `965fa34d9d051467d821e87c62806c8cea866e1478e3fc83f6bc5d20d2fe3f50`）で、隔離sourceとは異なります。隔離sourceとの差は393行追加、21行削除です。既存のMATC Node Cache patch（SHA-256 `53f97c44e82215fb3d230c0b90b024b126cd1dfde716a3ff3ec56c171c02c045`）は公式sourceにも現隔離sourceにも単独適用できず、異なる中間Phase24ソースを前提にしています。新しい単一パッチは隔離source全体の実体から生成しており、この中に現行cacheコードも含みます。

## 変更ファイルの機能分類

再構成パッチは公式5a8…に対する16ファイル差分です。ファイルごとのupstream blob、隔離ソースblob、raw SHA-256は `source_hash_manifest.json` を参照してください。

| 機能群 | ファイル | 差分の役割 |
|---|---|---|
| 熱solver / TES / assembly | `fem/src/modules/HeatSolve.F90` | Phase24のbulk/prism/tetra assembly、温度依存k、TES回路と計測、既定OFFのVector/Static/Dynamic K cache経路 |
| Vector/CSR assembly | `fem/src/CRSMatrix.F90`, `fem/src/ElementUtils.F90`, `fem/src/DefUtils.F90`, `fem/src/ListMatrix.F90` | CSR位置/glue経路と呼び出し補助、matrix utilityの連携 |
| Block / HYPRE / MPI solver | `fem/src/BlockSolve.F90`, `fem/src/SParIterSolver.F90`, `fem/src/SolverUtils.F90`, `fem/src/SolveHypre.c`, `fem/src/amgx.c`, `fem/src/ElmerSolver.F90` | mortar/block-Schur、HYPRE/MGR・MPI bridge、solver routingと実装群 |
| timestep / preconditioner | `fem/src/TimeIntegrate.F90`, `fem/src/VankaCreate.F90` | Phase24時間積分・前処理連携 |
| build / registry / metadata | `fem/src/CMakeLists.txt`, `fem/src/SOLVER.KEYWORDS`, `fem/src/ElementDescription.F90` | target登録、設定keyword、要素記述の連携 |

分類はファイル差分とソース内のPhase24/HYPRE/assembly呼び出しに基づくものです。各機能の個別性能・数値妥当性を今回新たに測ったものではありません。

## 再現パッチと適用順

1. Elmer sourceを公式commit `5a8de867068be0568f09af40fb90fee300dfbede` にcheckoutする。
2. `source_hash_manifest.json` の16項目すべてについてbase Git blobが一致することを確認する。
3. `git apply --check 01_phase24_source_from_release_26_1.patch` を実行する。
4. `git apply 01_phase24_source_from_release_26_1.patch` を実行する。
5. 16項目すべてのafter Git blobがmanifestと一致することを確認する。
6. `replay_source_and_build.ps1` のconfigure/build設定で隔離Release buildを作る。

歴史的 `elmer_native_source.patch` と旧MATC cache patchは再現順に含めません。前者は異なるHeatSolve preimageが必要で、後者は別の中間HeatSolve状態を必要とするためです。現在のsourceを丸ごと表す新パッチに対して、追加適用すると重複/競合します。

## Release build結果

新規build dir `C:\Users\Nao\AppData\Local\Temp\phase24-repro-release-build` でCMake configureとfull `cmake --build --parallel 6` が完了しました。全2,237工程が成功し、`HeatSolve.F90.obj`、`HeatSolve.dll`、`ElmerSolver_mpi.exe` が生成されました。生ログは同じフォルダの `configure.log` と `build_release.log` です。

構成はRelease `-O3`、MSYS2 GCC/GFortran 16.1.0、CMake 4.4.0、Ninja 1.13.2、OpenMP 5.2、MS-MPI 10.1.12498.52、HYPRE 2.33.0、OpenBLASです。依存実行ファイルとlibrariesのhashは `build_environment.json` にあります。

新規 `HeatSolve.dll` SHA-256は `9687e206fd81d15b6bdf2b45caabc7e5e26bb84d19fc7f6b6f67b02bbc392a9d`、新規solver SHA-256は `4c673dc0639f4cabc29e2ccc4b3dd43576f210a9bdced67a9be20761a50c3b6` です。旧既存Release DLLはそれぞれ `59384fbe...0643b` と `ec74448b...5778` で、byte一致しません。PE timestampも新旧で異なります（HeatSolve: 2026-10-12 00:43:42 対 2026-10-11 20:48:44）。タイムスタンプ以外を含むbinary差の全理由は今回分解していないため、bit-reproducible binaryとは主張しません。

## MATC AST安全性API設計

APIの設計案は `matc_ast_api_design_ja.md` に記録しました。MATC式評価処理は変更していません。単純文字列判定は使用せず、評価を実行しないparse-only AST validator、builtin純粋性許可表、ユーザー関数shadowing拒否、負例テストを分離した設計です。API実装と正例/負例runtime testは未実施です。

## Cache寿命設計

現行 `HeatSolver` のSAVE cacheと初期化/終了経路の監査、寿命設計案は `cache_lifecycle_design_ja.md` に記録しました。solver finalizerの実在とModel/Mesh generation hookは確定できていません。既定値OFFを維持し、今回のbuildでcacheを有効にしていません。

## 実施しなかった項目

- MATC AST API実装・負例runtime tests: **SKIP**
- cache lifecycle source修正・runtime再初期化/Model/Mesh差替えtests: **SKIP**
- DualTES/single-pixel 4-step比較、性能再測定: **SKIP**
- 75 ms計算: 実行なし
- 既存成功済み高速化証拠: 変更なし

今回確立したのはPhase24 **ソースとRelease buildの再現性**です。MATC cacheの汎用安全性や修正後の数値一致は確立していません。
