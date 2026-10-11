# 隔離Fortranビルド失敗の診断

## 観測

前回Ninja buildの失敗はFortran module依存情報の内容やMATCキャッシュの構文ではなく、GFortran frontendの起動環境で再現した。

```powershell
& 'G:\Elmer-Migration\msys64\ucrt64\lib\gcc\x86_64-w64-mingw32\16.1.0\f951.exe' --version
```

終了値は`-1073741511` (`0xC0000139`, `STATUS_ENTRYPOINT_NOT_FOUND`)。最小Fortran compileも同じDLL環境で終了コード1となった。Windows loaderが返した数値コードが観測値であり、最初の失敗時にGUI loaderの説明文やNinja stderr全文は保存されていない。

```powershell
& 'G:\Elmer-Migration\msys64\ucrt64\bin\gfortran.exe' -c 'C:\Users\Nao\AppData\Local\Temp\phase24-cache-runtime\tiny_fail.f90' -o 'tiny_fail.o'
```

## 設定の切り分け

対象ソースは`.F90`。Ninjaは`gfortran -cpp ... -E`による前処理後、`cmake -E cmake_ninja_depends`と`cmake_ninja_dyndep`でFortran依存を生成し、`-fpreprocessed -c`でcompileする。自由形式として処理される。ucrt64とMSYS2のDLL directoryをPATH先頭に設定した後、Vector Assembly OFFの対照buildとキャッシュ試作buildがDLL linkまで完了した。build logでは`-fopenmp`、Debug `-O0 -g3 -fbacktrace -fno-omit-frame-pointer`を確認した。

したがって再現した停止要因はFortran module生成やsource formではなくfrontendのDLL/ABI解決環境と判断する。完全な失敗時Ninja stderrが未保存である点は記録上の限界として残す。

## 証跡

- `build_failure_diagnosis.log`: frontend確認コマンドと終了値
- `build_off_control.log`: Vector Assembly OFF対照build
- `build_trial_debug.log`: 試作隔離Debug build
- 製品DLL、製品SIF、製品設定ファイルは変更していない
