# Source and build file locations

- `main.py` is the only Python file in the project root and the primary CLI.
- `src/` contains all maintained Python modules, analysis/preparation tools,
  support utilities, visualization commands, and PowerShell launchers.
- `src/fortran/` contains maintained Fortran user-function sources.
- `build/udf/` is for compiled Elmer user-function libraries (`.dll` / `.so`).
- `build/fortran/` is for compiler intermediates such as `.o` and `.mod` files.
- `tools/elmer-runtime/` is a local Elmer installation; it is excluded from Git.

Compile a user function from the project root after creating `build/udf/` if it
does not exist yet:

```powershell
New-Item -ItemType Directory -Force build/udf, build/fortran | Out-Null
elmerf90 src/fortran/tes_transient_heat_source.f90 -o build/udf/tes_transient_heat_source_t0.dll
```

`src/run.py` adds `build/udf/` to the native library lookup paths for Elmer runs.
