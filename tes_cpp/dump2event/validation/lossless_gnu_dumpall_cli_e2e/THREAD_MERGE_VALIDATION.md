# Threaded Dumpall merge validation status

A small 100-history run was launched with `$OMP=2` using the existing `phits_WinGfort_lossless_nohook.exe`. Inspection of the source makefile shows its default `USEOMP` is commented out; the `$OMP=2` input does not make that binary an OpenMP build. That run exited 0 and produced a valid single-thread-format Dumpall, but it is **not** evidence of threaded merge.

A fresh isolated copy of `C:\PHITS\src` was patched with the chunked cpcasdump writer and GNU-Fortran-16 compatibility rename, then built with:

```text
mingw32-make -f makefile ENVFLAGS=WinGfort USEOMP=true -j4
```

The build failed at `mod_ompparallel.o`; the assembler emitted `Error: unknown pseudo-op: .tls_common`. No OpenMP PHITS binary was produced. `C:\PHITS`, the original input, and the original Dumpall were not changed.

**Status:** multi-thread serialization merge order remains unvalidated. Do not infer it from the serial runs.
