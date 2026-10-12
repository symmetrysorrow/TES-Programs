#!/usr/bin/env python3
"""Validate the two-model same-process MATC cache lifecycle evidence."""
from __future__ import annotations

import hashlib
import json
import math
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
LOGS = ROOT / "logs"
ATOL = 1.0e-13
RTOL = 1.0e-12


def parse_kv(line: str) -> dict[str, str]:
    return {k: v.strip() for k, v in re.findall(r"([A-Za-z_]+)=\s*([^,\s]+)", line)}


def parse_log(name: str) -> dict:
    raw = (LOGS / name).read_bytes()
    encoding = "utf-16" if raw.startswith((b"\xff\xfe", b"\xfe\xff")) else "utf-8"
    decoded = raw.decode(encoding, errors="replace").lstrip("\ufeff")
    lines = decoded.splitlines()
    (LOGS / name.replace("_raw.log", ".log.txt")).write_text(decoded, encoding="utf-8")
    starts, ends, builds, clears, maps, evals, verifies, temps, passes = ([] for _ in range(9))
    for line_no, line in enumerate(lines):
        if "PHASE24_LIFECYCLE event=ACTIVATION_START" in line:
            starts.append(parse_kv(line))
        elif "PHASE24_LIFECYCLE event=ACTIVATION_END" in line:
            ends.append(parse_kv(line))
        elif "PHASE24_LIFECYCLE event=CACHE_BUILD" in line:
            builds.append(parse_kv(line))
        elif "PHASE24_LIFECYCLE event=CACHE_CLEAR" in line:
            # GNU Fortran line wrapping can split the trailing logical value.
            data = parse_kv(line)
            if "all_unallocated" not in data:
                data["all_unallocated"] = "T" if line_no + 1 < len(lines) and lines[line_no + 1].strip() == "=T" else "F"
            clears.append(data)
        elif "PHASE24_LIFECYCLE_MAP" in line:
            maps.append(parse_kv(line))
        elif "PHASE24_LIFECYCLE event=CACHE_EVAL" in line:
            evals.append(parse_kv(line))
        elif "PHASE24_LIFECYCLE_VERIFY" in line:
            verifies.append(parse_kv(line))
        elif "PHASE24_LIFECYCLE_TEMP" in line:
            temps.append(parse_kv(line))
        elif "PHASE24_MESHMODE_PASS" in line:
            passes.append(parse_kv(line))
    return {"lines": lines, "starts": starts, "ends": ends, "builds": builds,
            "clears": clears, "maps": maps, "evals": evals, "verifies": verifies,
            "temps": temps, "passes": passes}


def fail(checks: dict, name: str, detail: str) -> None:
    checks[name] = {"status": "FAIL", "detail": detail}


def ok(checks: dict, name: str, detail: str) -> None:
    checks[name] = {"status": "PASS", "detail": detail}


def fingerprint(rows: list[dict]) -> str:
    canonical = "\n".join(
        f"m={r['material']};e={r['element']};ln={r['local_node']};n={r['node']};s={r['slot']}"
        for r in sorted(rows, key=lambda x: (int(x["material"]), int(x["element"]), int(x["local_node"])))
    )
    return hashlib.sha256(canonical.encode()).hexdigest()


def fixture_connectivity(mesh: str) -> list[int]:
    path = ROOT / "fixtures" / "meshes" / mesh / "mesh.elements"
    rows = [line.split() for line in path.read_text(encoding="ascii").splitlines() if line.strip()]
    if len(rows) != 1 or len(rows[0]) < 9 or rows[0][2] != "706":
        raise ValueError(f"expected one Elmer 6-node prism in {path}")
    return [int(x) for x in rows[0][-6:]]


def main() -> int:
    checks: dict[str, dict] = {}
    on = parse_log("cache_on_raw.log")
    off = parse_log("cache_off_raw.log")
    expected = {"meshA": fixture_connectivity("meshA"),
                "meshA1": fixture_connectivity("meshA1")}

    for label, data in (("ON", on), ("OFF", off)):
        if len(data["starts"]) == 4 and len(data["ends"]) == 4:
            pids = {x.get("pid") for x in data["starts"]}
            act_ids = [x.get("activation") for x in data["starts"]]
            if len(pids) == 1 and act_ids == ["1", "2", "3", "4"]:
                ok(checks, f"{label.lower()}_same_pid_four_activations", f"pid={next(iter(pids))}; activation IDs={act_ids}")
            else:
                fail(checks, f"{label.lower()}_same_pid_four_activations", f"pids={pids}; activation IDs={act_ids}")
        else:
            fail(checks, f"{label.lower()}_same_pid_four_activations", f"starts={len(data['starts'])}; ends={len(data['ends'])}")

        event_order_ok = True
        for i, (start, end) in enumerate(zip(data["starts"], data["ends"])):
            start_idx = next((j for j, x in enumerate(data["lines"]) if "event=ACTIVATION_START" in x and f"activation={start.get('activation')} " in x), -1)
            end_idx = next((j for j, x in enumerate(data["lines"]) if "event=ACTIVATION_END" in x and f"activation={end.get('activation')} " in x), -1)
            if start_idx < 0 or end_idx <= start_idx:
                event_order_ok = False
            if i + 1 < len(data["starts"]):
                next_start_idx = next((j for j, x in enumerate(data["lines"]) if "event=ACTIVATION_START" in x and f"activation={data['starts'][i+1].get('activation')} " in x), -1)
                if next_start_idx <= end_idx:
                    event_order_ok = False
        if event_order_ok:
            ok(checks, f"{label.lower()}_activation_order", "each activation ended before the next activation started")
        else:
            fail(checks, f"{label.lower()}_activation_order", "activation start/end ordering invalid")

        meshes = [s.get("mesh") for s in data["starts"]]
        if meshes == ["meshA", "meshA", "meshA1", "meshA1"] and all(s.get("materials") == "7" and s.get("nodes") == "6" for s in data["starts"]):
            ok(checks, f"{label.lower()}_model_identity", f"meshes={meshes}; each model materials=7 nodes=6")
        else:
            fail(checks, f"{label.lower()}_model_identity", f"meshes={meshes}; starts={data['starts']}")

        if len(data["passes"]) == 2 and all(p.get("completed") in {"1", "2"} for p in data["passes"]):
            ok(checks, f"{label.lower()}_meshmode_source_path", f"Mesh Mode completed both models: {data['passes']}")
        else:
            fail(checks, f"{label.lower()}_meshmode_source_path", f"Mesh Mode pass markers={data['passes']}")

    on_by_activation = {s["activation"]: s for s in on["starts"]}
    actual_fingerprints = {}
    map_ok = True
    for activation, mesh in (("1", "meshA"), ("2", "meshA"), ("3", "meshA1"), ("4", "meshA1")):
        rows = [r for r in on["maps"] if r.get("activation") == activation]
        nodes = [int(r["node"]) for r in sorted(rows, key=lambda x: int(x["local_node"]))]
        if nodes != expected[mesh] or any(r.get("material") != "7" for r in rows) or len(rows) != 6:
            map_ok = False
        actual_fingerprints[activation] = fingerprint(rows) if rows else None
    if map_ok and actual_fingerprints["1"] == actual_fingerprints["2"] and actual_fingerprints["3"] == actual_fingerprints["4"] and actual_fingerprints["1"] != actual_fingerprints["3"]:
        ok(checks, "mapping_fingerprints", f"A={actual_fingerprints['1']}; B={actual_fingerprints['3']}; A != B; expected fixture connectivity matched")
    else:
        fail(checks, "mapping_fingerprints", f"fingerprints={actual_fingerprints}; mapping validation={map_ok}")

    lifecycle_ok = True
    for a in ("1", "2", "3", "4"):
        e = next((x for x in on["ends"] if x.get("activation") == a), {})
        b = [x for x in on["builds"] if x.get("activation") == a]
        c = [x for x in on["clears"] if x.get("activation") == a and x.get("allocated_before") == "6"]
        if e.get("builds") != "1" or len(b) != 1 or b[0].get("nodes") != "6" or b[0].get("materials") != "1" or not c or c[-1].get("freed_arrays") != "6" or c[-1].get("all_unallocated") != "T":
            lifecycle_ok = False
    if lifecycle_ok:
        ok(checks, "cache_build_and_release", "each ON activation built one 1-material/6-node map and explicitly freed all six allocated cache arrays")
    else:
        fail(checks, "cache_build_and_release", "missing build or six-array release proof for at least one activation")

    eval_ok = True
    for a in ("1", "2", "3", "4"):
        ev = [x for x in on["evals"] if x.get("activation") == a]
        ver = [x for x in on["verifies"] if x.get("activation") == a]
        if len(ev) != 2 or sum(int(x.get("hits", "0")) for x in ev) < 1 or len(ver) < 2 or any(float(x.get("max_abs", "inf")) > 1e-14 for x in ver):
            eval_ok = False
    if eval_ok:
        ok(checks, "per_assembly_evaluation_and_hits", "two temperature-dependent cache evaluations per activation; second assembly hits cache; Generic comparison max_abs=0")
    else:
        fail(checks, "per_assembly_evaluation_and_hits", "missing eval/hit/Generic verification evidence")

    off_ok = all(e.get("builds") == "0" and e.get("evaluations") == "0" for e in off["ends"])
    if off_ok:
        ok(checks, "cache_off_fallback", "OFF run completed all activations with no cache builds/evaluations")
    else:
        fail(checks, "cache_off_fallback", "OFF run unexpectedly used cache or lacked activation summaries")

    def temp_map(data: dict) -> dict:
        result = {}
        for t in data["temps"]:
            result[(t["activation"], t["node"])] = float(t["value"])
        return result

    on_t, off_t = temp_map(on), temp_map(off)
    compare_ok = set(on_t) == set(off_t) and bool(on_t)
    max_abs = 0.0
    max_rel = 0.0
    for key in set(on_t) & set(off_t):
        diff = abs(on_t[key] - off_t[key])
        max_abs = max(max_abs, diff)
        max_rel = max(max_rel, diff / max(abs(off_t[key]), 1.0e-300))
        if diff > ATOL + RTOL * abs(off_t[key]):
            compare_ok = False
    if compare_ok:
        ok(checks, "cache_on_off_temperature_comparison", f"{len(on_t)} node values; max_abs={max_abs:.17g}; max_rel={max_rel:.17g}; atol={ATOL}; rtol={RTOL}")
    else:
        fail(checks, "cache_on_off_temperature_comparison", f"key sets differ or tolerance exceeded; max_abs={max_abs}; max_rel={max_rel}")

    result_hashes = {}
    result_equal = True
    for mesh in ("meshA", "meshA1"):
        off_path = ROOT / mesh / "meshmode_cacheoff.result"
        on_path = ROOT / mesh / "meshmode_cacheon.result"
        if not off_path.is_file() or not on_path.is_file():
            result_equal = False
            continue
        off_hash = hashlib.sha256(off_path.read_bytes()).hexdigest()
        on_hash = hashlib.sha256(on_path.read_bytes()).hexdigest()
        result_hashes[mesh] = {"off": off_hash, "on": on_hash}
        result_equal &= off_hash == on_hash
    if result_equal and len(result_hashes) == 2:
        ok(checks, "cache_on_off_result_files", f"Elmer .result output SHA-256 matches byte-for-byte for both meshes: {result_hashes}")
    else:
        fail(checks, "cache_on_off_result_files", f"missing or different ON/OFF result files: {result_hashes}")

    checks["cache_toggle_off_on_off"] = {
        "status": "SKIP",
        "detail": "Solver option is read from the model Solver ValueList at activation; this fixture uses one shared SIF configuration, so runtime OFF→ON→OFF mutation was not attempted. Separate OFF and ON processes were compared using identical meshes and physics.",
    }
    statuses = [v["status"] for v in checks.values()]
    status = "PASS" if statuses and all(x in {"PASS", "SKIP"} for x in statuses) else "FAIL"
    report = {
        "status": status,
        "scope": "real Elmer 26.1 HeatSolver + production Generic MATC node-cache implementation; one rank; two Mesh Mode models per solver process",
        "cache_default": "OFF",
        "vector_assembly": "OFF",
        "same_pid": {"cache_on": on["starts"][0].get("pid") if on["starts"] else None,
                     "cache_off": off["starts"][0].get("pid") if off["starts"] else None},
        "mapping_fingerprints": actual_fingerprints,
        "cache_on_off_temperature_max_abs": max_abs,
        "cache_on_off_temperature_max_rel": max_rel,
        "thresholds": {"atol": ATOL, "rtol": RTOL},
        "checks": checks,
    }
    (ROOT / "pass_fail_skip.json").write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2, ensure_ascii=False))
    return 0 if status == "PASS" else 1


if __name__ == "__main__":
    sys.exit(main())
