from pathlib import Path
import hashlib, json, subprocess

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[5]
SOURCE = Path(r"C:\Users\Nao\AppData\Local\Temp\phase24-source-replay-worktree")
DIFF = Path(r"C:\Users\Nao\AppData\Local\Temp\phase24-02-diff")
BUILD = Path(r"C:\Users\Nao\AppData\Local\Temp\phase24-repro-release-build")
UPSTREAM = Path(r"C:\Users\Nao\AppData\Local\Temp\elmer-release-26.1-full")
PATCH01 = HERE.parent / "phase24-source-repro" / "01_phase24_source_from_release_26_1.patch"
PATCH02 = HERE / "02_matc_ast_purity_validator.patch"
MANIFEST01 = json.loads((PATCH01.parent / "source_hash_manifest.json").read_text())

def sha(path):
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()

def git(*args):
    return subprocess.check_output(["git", *map(str, args)], text=True).strip()

changed = ["fem/src/modules/HeatSolve.F90", "matc/src/elmer/fnames.h",
          "matc/src/matc.c", "matc/src/parser.c"]
patch_files = []
for path in changed:
    before = git("-C", DIFF, "rev-parse", f":{path}")
    after = git("-C", SOURCE, "hash-object", f"--path={path}", SOURCE / path)
    patch_files.append({"path": path, "before_git_blob": before, "after_git_blob": after,
                        "before_raw_sha256": sha(UPSTREAM / path) if path != "fem/src/modules/HeatSolve.F90" else
                            next(x["target_raw_sha256"] for x in MANIFEST01["files"] if x["path"] == path),
                        "after_raw_sha256": sha(SOURCE / path)})

source_replay = []
for f in MANIFEST01["files"]:
    actual = git("-C", DIFF, "rev-parse", f":{f['path']}")
    source_replay.append({"path": f["path"], "target_git_blob": f["target_git_blob"],
                          "replayed_01_index_git_blob": actual,
                          "matches_01": actual == f["target_git_blob"]})
artifacts = {
    "ElmerSolver_mpi.exe": BUILD / "fem/src/ElmerSolver_mpi.exe",
    "HeatSolve.dll": BUILD / "fem/src/modules/HeatSolve.dll",
    "libmatc.dll": BUILD / "matc/src/libmatc.dll",
    "libHYPRE.dll": Path(r"G:\Elmer-Migration\msys64\ucrt64\bin\libHYPRE.dll"),
    "msmpi.dll": Path(r"C:\Windows\System32\msmpi.dll"),
}
binary_hashes = {name: {"path": str(path), "sha256": sha(path)} for name, path in artifacts.items() if path.is_file()}
cache = (BUILD / "CMakeCache.txt").read_text(errors="replace").splitlines()
keys = ("CMAKE_BUILD_TYPE:STRING", "CMAKE_C_COMPILER:STRING", "CMAKE_CXX_COMPILER:FILEPATH",
        "CMAKE_Fortran_COMPILER:FILEPATH", "CMAKE_MAKE_PROGRAM:FILEPATH", "WITH_MPI:BOOL",
        "WITH_Hypre:BOOL", "WITH_OpenMP:BOOL")
cmake = {k.split(":")[0]: next((line.split("=",1)[1] for line in cache if line.startswith(k+"=")), None) for k in keys}
result = {
    "schema": "phase24-matc-ast-source-build-v1",
    "upstream": {"repository": "ElmerCSC/elmerfem", "commit": "5a8de867068be0568f09af40fb90fee300dfbede"},
    "patches": [
        {"path": str(PATCH01.relative_to(REPO)).replace("\\", "/"), "sha256": sha(PATCH01),
         "order": 1, "matches_existing_manifest": sha(PATCH01) == MANIFEST01["patch"]["sha256"]},
        {"path": str(PATCH02.relative_to(REPO)).replace("\\", "/"), "sha256": sha(PATCH02),
         "order": 2, "file_count": len(patch_files)},
    ],
    "phase24_01_files": source_replay,
    "phase24_01_all_match": all(x["matches_01"] for x in source_replay),
    "patch_02_files": patch_files,
    "build": {"configuration": cmake, "compiler": {
        "gcc": subprocess.check_output([r"G:\Elmer-Migration\msys64\ucrt64\bin\gcc.exe", "--version"], text=True).splitlines()[0],
        "gfortran": subprocess.check_output([r"G:\Elmer-Migration\msys64\ucrt64\bin\gfortran.exe", "--version"], text=True).splitlines()[0],
        "cmake": subprocess.check_output(["cmake", "--version"], text=True).splitlines()[0],
        "ninja": subprocess.check_output(["ninja", "--version"], text=True).strip()},
        "clean_release_build_log": "build_release_ast.log",
        "final_incremental_matc_build_log": "build_matc_ast_final.log",
        "binaries": binary_hashes},
    "tests": {"ast_unit_test_exit": 0, "ast_unit_test_log": "matc_ast_purity_test_final.log",
              "final_four_step_runs": "runtime_runs.json", "dualtes_matrix_rhs": "matrix_rhs_comparison.json"},
    "runtime_summary": {
        "vector_assembly": False,
        "runs": [{k: r[k] for k in ("case", "mode", "returncode", "wall_seconds", "ranks", "steps", "cache",
                                     "solver_sha256", "heat_solve_sha256")}
                 for r in json.loads((HERE / "runtime_runs.json").read_text())["records"]],
        "pre_ast_comparison": "pre_ast_runtime_comparison.json",
        "matrix_rhs_comparison": json.loads((HERE / "matrix_rhs_comparison.json").read_text()),
        "source_replay_log": "source_replay_apply.log"},
    "artifact_inventory": "artifact_inventory.json",
}
(HERE / "02_source_hash_manifest.json").write_text(json.dumps({
    "schema": "phase24-matc-ast-patch-source-hashes-v1",
    "patch": {"path": "02_matc_ast_purity_validator.patch", "sha256": sha(PATCH02), "order": 2},
    "files": patch_files,
    "binaries": binary_hashes}, indent=2) + "\n", encoding="utf-8")
(HERE / "evidence_manifest.json").write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
print(json.dumps({"phase24_01_all_match": result["phase24_01_all_match"],
                  "patch_sha256": sha(PATCH02), "changed_files": patch_files,
                  "binary_hashes": binary_hashes}, indent=2))
