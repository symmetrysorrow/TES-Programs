#!/usr/bin/env python3
"""Run four-step cache OFF/ON cases against the reconstructed Phase24 Release."""
from __future__ import annotations
import argparse, hashlib, json, os, re, shutil, subprocess, time
from pathlib import Path

HERE = Path(__file__).resolve().parent
OLD_EVIDENCE = HERE.parent / "matc-cache-release"
PARENT = HERE.parent / "bulk-matc-cache"
RUNTIME = Path(r"C:\Users\Nao\AppData\Local\Temp\phase24-cache-runtime")
BUILD = Path(r"C:\Users\Nao\AppData\Local\Temp\phase24-repro-release-build")
SOLVER = BUILD / "fem" / "src" / "ElmerSolver_mpi.exe"
MPI = Path(r"C:\Program Files\Microsoft MPI\Bin\mpiexec.exe")
MESHES = RUNTIME / "workspace" / "work" / "meshes"
CASES = [
    ("dualtes", "mesh_toml_95333801e1a3", "dualtes_off_1.sif", "dualtes_on_1.sif", [PARENT / f"phase24_short_ab_dualtes_vector_{i}.state" for i in (1, 2)]),
    ("singlepixel", "mesh_hybrid_fullconf_h8", "singlepixel_release_off_6rank.sif", "singlepixel_release_on_6rank.sif", [HERE.parent / "phase24_short_ab_singlepixel_vector_1.state"]),
]
PATHS = [RUNTIME / "build" / "udf", BUILD / "fem" / "src" / "modules", BUILD / "fem" / "src",
         BUILD / "matc" / "src", BUILD / "fhutiter" / "src", BUILD / "fem" / "src" / "binio",
         BUILD / "mathlibs" / "src" / "arpack", BUILD / "mathlibs" / "src" / "parpack",
         RUNTIME / "bin", RUNTIME / "share" / "elmersolver" / "lib",
         Path(r"G:\Elmer-Migration\msys64\ucrt64\bin"), Path(r"G:\Elmer-Migration\msys64\usr\bin"),
         Path(r"G:\Elmer-Migration\elmer\bin"), MPI.parent]

def sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()

def parse_result(path: Path):
    lines = path.read_text(encoding="ascii", errors="replace").splitlines()
    result = {}; i = 0
    while i < len(lines):
        if not lines[i].startswith("Time:"):
            i += 1; continue
        step = int(lines[i].split()[1]); i += 1
        variable = lines[i].strip(); i += 1
        if not lines[i].startswith("Perm:"):
            raise ValueError(f"{path}: expected Perm at line {i+1}")
        count = int(lines[i].split()[2]); i += 1
        ids = [int(lines[i+j].split()[1]) for j in range(count)]; i += count
        values = [float(lines[i+j].strip()) for j in range(count)]; i += count
        result[step] = (variable, dict(zip(ids, values)))
    return result

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--case", choices=[x[0] for x in CASES])
    parser.add_argument("--mode", choices=["off", "on"])
    args = parser.parse_args()
    if not SOLVER.is_file() or not MPI.is_file():
        raise SystemExit("reconstructed Release solver or MS-MPI is missing")
    env = os.environ.copy()
    env.update(ELMER_HOME=str(RUNTIME), OMP_NUM_THREADS="1", OMP_DYNAMIC="FALSE", PHASE24_HYPRE_REPORT="1")
    env["PATH"] = os.pathsep.join(map(str, PATHS)) + os.pathsep + env.get("PATH", "")
    previous = HERE / "runtime_partial_runs.json"
    records = json.loads(previous.read_text(encoding="utf-8")).get("records", []) if (args.case or args.mode) and previous.exists() else []
    selected_cases = [c for c in CASES if not args.case or c[0] == args.case]
    if args.case or args.mode:
        records = [r for r in records if not ((not args.case or r["case"] == args.case) and (not args.mode or r["mode"] == args.mode))]
    for case, mesh_name, off_sif, on_sif, states in selected_cases:
        mesh = MESHES / mesh_name
        selected_modes = [(m, s) for m, s in (("off", off_sif), ("on", on_sif)) if not args.mode or m == args.mode]
        for mode, input_sif in selected_modes:
            for state in states:
                if not state.is_file():
                    raise SystemExit(f"missing restart state: {state}")
                shutil.copy2(state, mesh / state.name)
            raw = (OLD_EVIDENCE / input_sif).read_text(encoding="utf-8-sig")
            prefix = f"ast_{case}_{mode}"
            # Give every solver/data output a unique name while preserving the
            # mesh, restart, four-step timing and all physical settings.
            raw = re.sub(r"(?m)^\s*Output File\s*=.*$", f"  Output File = {prefix}.result", raw, count=1)
            if case == "dualtes":
                raw = raw.replace(f"dualtes_release_{mode}_1", prefix)
                for tes in (1, 2):
                    raw = re.sub(rf'(?m)^(\s*TES\s+{tes}\s+Series File\s*=\s*String\s+")[^"\r\n]*(".*)$',
                                 rf'\g<1>{prefix}_{tes}_series.csv\2', raw)
                    raw = re.sub(rf'(?m)^(\s*TES\s+{tes}\s+Iteration Series File\s*=\s*String\s+")[^"\r\n]*(".*)$',
                                 rf'\g<1>{prefix}_{tes}_iterations.csv\2', raw)
            else:
                raw = raw.replace(f"singlepixel_release_{mode}", prefix)
                raw = re.sub(r'(?m)^(\s*TES Series File\s*=\s*String\s+")[^"\r\n]*(".*)$',
                             rf'\g<1>{prefix}_1_series.csv\2', raw)
                raw = re.sub(r'(?m)^(\s*TES Iteration Series File\s*=\s*String\s+")[^"\r\n]*(".*)$',
                             rf'\g<1>{prefix}_1_iterations.csv\2', raw)
            raw = re.sub(r"(?m)^\s*Phase24 Outer HeatSolve Capture Prefix\s*=.*$", f'  Phase24 Outer HeatSolve Capture Prefix = String "{prefix}"', raw)
            sif = HERE / f"{prefix}.sif"
            sif.write_text(raw, encoding="utf-8", newline="\n")
            shutil.copy2(sif, RUNTIME / sif.name)
            log = HERE / f"{prefix}.solver.log"
            t0 = time.perf_counter()
            with log.open("w", encoding="utf-8", newline="\n") as out:
                proc = subprocess.run([str(MPI), "-n", "6", str(SOLVER), sif.name], cwd=RUNTIME, env=env,
                                      stdout=out, stderr=subprocess.STDOUT, check=False)
            wall = time.perf_counter() - t0
            text = log.read_text(encoding="utf-8", errors="replace")
            outputs = []
            for tes in (range(1, 3) if case == "dualtes" else range(1, 2)):
                for kind in ("series", "iterations"):
                    p = RUNTIME / f"{prefix}_{tes}_{kind}.csv"
                    if p.exists():
                        dest = HERE / p.name; shutil.copy2(p, dest)
                        outputs.append({"file": dest.name, "sha256": sha(dest)})
            norms = [float(x) for x in re.findall(r"HeatSolve:\s+Result Norm\s*:\s*([0-9.Ee+-]+)", text)]
            pcg = [[int(i), float(r)] for i, r in re.findall(r"Phase24 PCG iterations=(\d+) final_rel_res=([0-9.Ee+-]+)", text)]
            result_base = RUNTIME / mesh_name / f"{prefix}.result"
            result_shards = [p for p in [result_base, *sorted((RUNTIME / mesh_name).glob(result_base.name + ".*"))] if p.is_file()]
            outputs += [{"file": p.name, "bytes": p.stat().st_size, "sha256": sha(p)} for p in result_shards]
            rec = {"case": case, "mode": mode, "returncode": proc.returncode, "wall_seconds": wall,
                   "ranks": 6, "steps": 4, "vector_assembly": False, "cache": mode == "on",
                   "solver_sha256": sha(SOLVER), "heat_solve_sha256": sha(BUILD / "fem" / "src" / "modules" / "HeatSolve.dll"),
                   "sif_sha256": sha(sif), "log": log.name, "result_norms": norms,
                   "hypre_pcg": pcg, "outputs": outputs}
            records.append(rec)
            (HERE / "runtime_partial_runs.json").write_text(json.dumps({"build": "reconstructed official-source Phase24 Release", "records": records}, indent=2)+"\n", encoding="utf-8")
            print(json.dumps(rec, separators=(",", ":")), flush=True)
            if proc.returncode != 0:
                raise SystemExit(f"solver failed: {case} {mode} rc={proc.returncode}")
    # Full nodal temperature comparisons for each reconstructed OFF/ON pair.
    comparisons = {}
    for case, mesh_name, *_ in CASES:
        mesh = RUNTIME / mesh_name; per = []; hashes = []
        for rank in range(6):
            a = mesh / f"ast_{case}_off.result.{rank}"; b = mesh / f"ast_{case}_on.result.{rank}"
            aa, bb = parse_result(a), parse_result(b)
            hashes.extend([{"file": p.name, "bytes": p.stat().st_size, "sha256": sha(p)} for p in (a, b)])
            if aa.keys() != bb.keys(): raise SystemExit(f"step mismatch: {case}, rank {rank}")
            for step in sorted(aa):
                av, bv = aa[step][1], bb[step][1]
                if av.keys() != bv.keys(): raise SystemExit(f"node mismatch: {case}, rank {rank}, step {step}")
                triples = [(av[k], bv[k], abs(av[k]-bv[k])) for k in av]
                per.append({"rank": rank, "timestep": step, "nodes": len(av),
                            "max_abs_K": max((x[2] for x in triples), default=0.0),
                            "max_rel": max((x[2]/max(abs(x[0]),abs(x[1]),1e-300) for x in triples), default=0.0),
                            "different_values": sum(x[2] != 0.0 for x in triples)})
        comparisons[case] = {"ranks": 6, "steps": 4, "temperature": "all nodes, rank shards", "per_rank_timestep": per,
                             "max_abs_K": max(x["max_abs_K"] for x in per), "max_rel": max(x["max_rel"] for x in per),
                             "total_different_values": sum(x["different_values"] for x in per), "result_shard_hashes": hashes}
    (HERE / "runtime_runs.json").write_text(json.dumps({"build": "reconstructed official-source Phase24 Release",
       "records": records, "temperature_comparisons": comparisons}, indent=2)+"\n", encoding="utf-8")
    return 0

if __name__ == "__main__": raise SystemExit(main())
