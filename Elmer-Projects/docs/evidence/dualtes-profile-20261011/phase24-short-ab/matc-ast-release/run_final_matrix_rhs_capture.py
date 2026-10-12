from pathlib import Path
import hashlib, json, os, re, shutil, subprocess, time

HERE = Path(__file__).resolve().parent
SRC = HERE.parent / "matc-cache-release"
PARENT = HERE.parent / "bulk-matc-cache"
RUNTIME = Path(r"C:\Users\Nao\AppData\Local\Temp\phase24-cache-runtime")
BUILD = Path(r"C:\Users\Nao\AppData\Local\Temp\phase24-repro-release-build")
MESH = RUNTIME / "workspace" / "work" / "meshes" / "mesh_toml_95333801e1a3"
MPI = Path(r"C:\Program Files\Microsoft MPI\Bin\mpiexec.exe")
SOLVER = BUILD / "fem" / "src" / "ElmerSolver_mpi.exe"
STATES = [PARENT / f"phase24_short_ab_dualtes_vector_{i}.state" for i in (1, 2)]
PATHS = [RUNTIME / "build" / "udf", BUILD / "fem" / "src" / "modules", BUILD / "fem" / "src",
         BUILD / "matc" / "src", BUILD / "fhutiter" / "src", BUILD / "fem" / "src" / "binio",
         BUILD / "mathlibs" / "src" / "arpack", BUILD / "mathlibs" / "src" / "parpack",
         RUNTIME / "bin", RUNTIME / "share" / "elmersolver" / "lib",
         Path(r"G:\Elmer-Migration\msys64\ucrt64\bin"), Path(r"G:\Elmer-Migration\msys64\usr\bin"),
         Path(r"G:\Elmer-Migration\elmer\bin"), MPI.parent]
env = os.environ.copy()
env.update(ELMER_HOME=str(RUNTIME), OMP_NUM_THREADS="1", OMP_DYNAMIC="FALSE", PHASE24_HYPRE_REPORT="1")
env["PATH"] = os.pathsep.join(map(str, PATHS)) + os.pathsep + env.get("PATH", "")

def sha(path):
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()

runs = []
for mode in ("off", "on"):
    label = f"final_matcap_{mode}"
    sif_text = (SRC / f"{mode}_ts1.sif").read_text(encoding="utf-8-sig")
    sif_text = sif_text.replace("Phase24 Outer HeatSolve Capture = Logical False",
                                "Phase24 Outer HeatSolve Capture = Logical True")
    sif_text = re.sub(r"(?m)^\s*Phase24 Outer HeatSolve Capture Timestep\s*=.*$",
                      "  Phase24 Outer HeatSolve Capture Timestep = Integer 1", sif_text)
    sif_text = re.sub(r"(?m)^\s*Phase24 Outer HeatSolve Capture Max Iterations\s*=.*$",
                      "  Phase24 Outer HeatSolve Capture Max Iterations = Integer 120", sif_text)
    sif_text = re.sub(r"(?m)^\s*Phase24 Outer HeatSolve Capture Prefix\s*=.*$",
                      f'  Phase24 Outer HeatSolve Capture Prefix = String "{label}"', sif_text)
    sif_text = re.sub(r"(?m)^\s*Output File\s*=.*$", f"  Output File = {label}.result", sif_text, count=1)
    sif_text = re.sub(r'(?m)^(\s*TES\s+\d+\s+Series File\s*=\s*String\s+")[^"\r\n]*(".*)$',
                      rf"\g<1>{label}_series\g<2>", sif_text)
    sif_text = re.sub(r'(?m)^(\s*TES\s+\d+\s+Iteration Series File\s*=\s*String\s+")[^"\r\n]*(".*)$',
                      rf"\g<1>{label}_iterations\g<2>", sif_text)
    sif = HERE / f"{label}.sif"
    sif.write_text(sif_text, encoding="utf-8", newline="\n")
    shutil.copy2(sif, RUNTIME / sif.name)
    capture_dir = RUNTIME / label
    for nl in range(1, 121):
        (capture_dir / f"ts0001_nl{nl:04d}").mkdir(parents=True, exist_ok=True)
    for state in STATES:
        shutil.copy2(state, MESH / state.name)
    log = HERE / f"{label}.solver.log"
    start = time.perf_counter()
    with log.open("w", encoding="utf-8", newline="\n") as stream:
        result = subprocess.run([str(MPI), "-n", "6", str(SOLVER), sif.name], cwd=RUNTIME,
                                env=env, stdout=stream, stderr=subprocess.STDOUT)
    files = sorted(p for p in capture_dir.rglob("*") if p.is_file())
    record = {"label": label, "mode": mode, "timestep": 1, "returncode": result.returncode,
              "wall_seconds": time.perf_counter() - start, "solver_sha256": sha(SOLVER),
              "heat_solve_sha256": sha(BUILD / "fem" / "src" / "modules" / "HeatSolve.dll"),
              "log": log.name, "capture_files": [
                  {"relative": p.relative_to(capture_dir).as_posix(), "bytes": p.stat().st_size, "sha256": sha(p)}
                  for p in files]}
    runs.append(record)
    (HERE / "matrix_rhs_capture_runs.json").write_text(json.dumps({
        "build": "reconstructed Phase24 Release with AST purity patch", "ranks": 6,
        "vector_assembly": False,
        "capture": "DualTES timestep 1, all nonlinear assemblies before/after DefaultSolve",
        "runs": runs}, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(record), flush=True)
    if result.returncode:
        raise SystemExit(f"{label} failed: {result.returncode}")

