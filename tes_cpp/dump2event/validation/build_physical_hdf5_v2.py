#!/usr/bin/env python3
"""Stream the validated Pb-only experimental ledger and raw tracks into HDF5 v2.

This is an opt-in research utility. It does not call the ordinary converter and
does not claim support for physical processes outside the existing Pb ledger.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import itertools
import json
import re
import subprocess
import sys
import time
import tempfile
from pathlib import Path

# Permit the installed/standalone validation driver to import the repository's
# small HDF5 v2 writer without requiring callers to edit PYTHONPATH.
_python_root = Path(__file__).resolve().parents[2] / "python"
if _python_root.is_dir() and str(_python_root) not in sys.path:
    sys.path.insert(0, str(_python_root))

from tes_cpp.event_hdf5_v2 import EventV2Writer, iter_deposits, iter_events


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def canonical_json(value: object) -> bytes:
    """Stable UTF-8 representation for semantic configuration hashes."""
    return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=False).encode("utf-8")


def source_definition(source_lines: list[str]) -> dict:
    """Parse only the explicitly supported unit-intensity source definitions."""
    values: dict[str, str] = {}
    energies: list[tuple[float, float]] = []
    for line in source_lines:
        if "=" in line:
            key, value = line.split("=", 1)
            values[re.sub(r"\s+", "", key).lower()] = value.strip()
            continue
        fields = line.split()
        if len(fields) == 2:
            try:
                energies.append((float(fields[0]), float(fields[1])))
            except ValueError:
                continue
    allowed = (([1.173, 1.332], "co60-lines"), ([0.511], "mono-gamma-511keV"),
               ([2.0], "mono-gamma-2MeV"))
    match = None
    for expected, mode in allowed:
        if len(energies) == len(expected) and all(
                abs(got[0] - want) < 1e-12 and abs(got[1] - 1.0) < 1e-12
                for got, want in zip(energies, expected)):
            match = (expected, mode)
            break
    if match is None:
        raise ValueError("source must be unit-intensity Co-60 lines or mono gamma at 0.511/2.0 MeV")
    expected, mode = match
    if values.get("totfact", "").split()[:1] not in (["1"], ["1.0"], ["1.00"]):
        raise ValueError("source must set unit totfact=1")
    if mode != "co60-lines":
        if values.get("proj", "").lower() != "photon":
            raise ValueError("mono gamma source must set proj=photon")
        if values.get("e-type", "") != "8" or values.get("ne", "") != "1":
            raise ValueError("mono gamma source must set e-type=8 and ne=1")
        if values.get("wgt", "1").split()[:1] not in (["1"], ["1.0"], ["1.00"]):
            raise ValueError("mono gamma source must use unit wgt=1")
    return {"source_mode": mode, "source_particle": "photon", "source_energies_MeV": expected}


def resolved_physics(settings: dict, source_input: Path) -> dict:
    """Return canonical effective model settings, excluding paths and byte hashes."""
    text = source_input.read_text(encoding="utf-8-sig", errors="replace")
    sections: dict[str, list[str]] = {}
    active = ""
    for raw in text.splitlines():
        line = raw.split("#", 1)[0].strip()
        if not line:
            continue
        if line.startswith("[") and line.endswith("]"):
            active = re.sub(r"\s+", "", line[1:-1]).lower()
            sections.setdefault(active, [])
        elif active:
            sections[active].append(line)
    return {
        "phits_version": str(settings["phits_version"]),
        "region_id": int(settings["region_id"]),
        "region_material": str(settings["region_material"]).upper(),
        "dumpall_material_id": int(settings["dumpall_material_id"]),
        "egs5_medium": int(settings["egs5_medium"]),
        "electron_cutoff_MeV": float(settings["electron_cutoff_MeV"]),
        "positron_cutoff_MeV": float(settings["positron_cutoff_MeV"]),
        "weight_mode": str(settings["weight_mode"]).lower(),
        "energy_unit": str(settings["energy_unit"]),
        "deposit_mode": str(settings["deposit_mode"]),
        "algorithm_version": "pb-egs5-exp-v2-positron",
        **source_definition(sections.get("source", [])),
    }


def physics_config(source_input: Path, config_path: Path | None = None):
    source_text = source_input.read_text(encoding="utf-8-sig", errors="replace")
    seed_match = re.search(r"(?im)^\s*bitrseed\s*=\s*([01]{64})\s*(?:#.*)?$", source_text)
    if not seed_match:
        raise ValueError("PHITS source input must declare a 64-bit binary bitrseed for reproducible provenance")
    sections: dict[str, list[str]] = {}
    active = ""
    for raw in source_text.splitlines():
        line = raw.split("#", 1)[0].strip()
        if not line:
            continue
        if line.startswith("[") and line.endswith("]"):
            active = re.sub(r"\s+", "", line[1:-1]).lower()
            sections.setdefault(active, [])
        elif active:
            sections[active].append(line)
    source = source_definition(sections.get("source", []))
    return {
        "phits_version": "3.341", "region_id": 1000, "region_material": "PB",
        "dumpall_material_id": 2, "egs5_medium": 2, "electron_cutoff_MeV": 0.001,
        "weights": "unit-weight validated slice only", "diagnostic_mode": True,
        **source,
        "unsupported_processes": "stored UNRESOLVED; never old-formula fallback",
        "source_input_sha256": sha256(source_input),
        "source_seed_bits": seed_match.group(1),
        "position_semantics": "raw transport chords and raw child point; scorer representative coordinate kept separate",
    }


def provenance_metadata(physics_settings: dict | None, config_path: Path | None,
                        source_input: Path, dumpall: Path, stopping_table: Path,
                        converter_executable: Path | None = None,
                        reconstructor: Path | None = None,
                        track_exporter: Path | None = None) -> dict:
    """Build distinct byte-level and semantic hashes for the HDF5 v2 provenance."""
    settings = physics_settings or {
        "phits_version": "3.341", "region_id": 1000, "region_material": "Pb",
        "dumpall_material_id": 2, "egs5_medium": 2, "electron_cutoff_MeV": 0.001,
        "positron_cutoff_MeV": 0.001, "weight_mode": "unit", "energy_unit": "MeV",
        "deposit_mode": "t-deposit-output-deposit",
    }
    resolved = resolved_physics(settings, source_input)
    resolved_hash = hashlib.sha256(canonical_json(resolved)).hexdigest()
    input_hash = sha256(source_input)
    physics = physics_config(source_input)
    metadata = {
        "provenance_schema_revision": 2,
        # Kept as a v2 compatibility alias, now defined identically to resolved_physics_sha256.
        "physics_config_sha256": resolved_hash,
        "resolved_physics_sha256": resolved_hash,
        "phits_input_sha256": input_hash,
        "source_input_sha256": input_hash,
        "dumpall_sha256": sha256(dumpall),
        "stopping_table_sha256": sha256(stopping_table),
        "phits_version": str(settings["phits_version"]),
        "resolved_physics_json": canonical_json(resolved).decode("utf-8"),
        "source_mode": physics["source_mode"],
        "source_energy_MeV": ",".join(format(x, ".12g") for x in physics["source_energies_MeV"]),
        "reconstruction_algorithm_version": "pb-egs5-exp-v2-positron",
        "validation_scope": "PHITS-3.341; Pb region 1000; unit-weight; validated EGS5 process subset",
        "tally_semantics": "PHITS [t-deposit] unit=3 output=deposit; transport-step position remains a segment",
        "driver_sha256": sha256(Path(__file__).resolve()),
    }
    if config_path is not None:
        metadata["physics_contract_sha256"] = sha256(config_path)
        contract = json.loads(config_path.read_text(encoding="utf-8"))
        phits_executable_sha = contract.get("phits_executable_sha256")
        if phits_executable_sha is not None:
            if not re.fullmatch(r"[0-9a-fA-F]{64}", str(phits_executable_sha)):
                raise ValueError("phits_executable_sha256 in physics contract must be a 64-character SHA-256 digest")
            metadata["phits_executable_sha256"] = str(phits_executable_sha).lower()
    for key, path in (("converter_executable_sha256", converter_executable),
                      ("reconstructor_executable_sha256", reconstructor),
                      ("transport_exporter_sha256", track_exporter)):
        if path is not None:
            metadata[key] = sha256(path)
    return metadata


def validate_source_input(path: Path) -> None:
    """Fail closed unless the input explicitly matches the validated Pb slice."""
    text = path.read_text(encoding="utf-8-sig", errors="replace")
    active = ""
    sections: dict[str, list[str]] = {}
    for raw in text.splitlines():
        line = raw.split("#", 1)[0].strip()
        if not line:
            continue
        if line.startswith("[") and line.endswith("]"):
            active = re.sub(r"\s+", "", line[1:-1]).lower()
            sections.setdefault(active, [])
        elif active:
            sections[active].append(line.lower())
    def values(section: str) -> dict[str, str]:
        found = {}
        for line in sections.get(section, []):
            if "=" not in line:
                continue
            key, value = line.split("=", 1)
            found[re.sub(r"\s+", "", key)] = value.strip()
        return found
    params = values("parameters")
    required_params = {"icntl": "0", "maxbch": "1",
                       "negs": "1", "e-mode": "2", "emin(12)": "0.001",
                       "emin(13)": "0.001", "emin(14)": "0.001"}
    for key, expected in required_params.items():
        got = params.get(key)
        if got is None or float(got.split()[0]) != float(expected):
            raise ValueError(f"source input must set [{key}]={expected}; found {got!r}")
    seed = params.get("bitrseed", "").split()[0] if params.get("bitrseed") else ""
    if not re.fullmatch(r"[01]{64}", seed):
        raise ValueError("source input must set an explicit 64-bit binary bitrseed")
    material_lines = sections.get("material", [])
    if not any(re.match(r"mat\s*\[\s*2\s*\]\s+pb\s+1(?:\s|$)", line) for line in material_lines):
        raise ValueError("source input must define mat[2] as pure Pb")
    cell_lines = sections.get("cell", [])
    if not any(re.match(r"1000\s+2\s+-9\.78(?:\s|$)", line) for line in cell_lines):
        raise ValueError("source input must define region 1000 with Dumpall material 2 and Pb density -9.78")
    source_definition(sections.get("source", []))
    tally = values("t-deposit")
    expected_tally = {"mesh": "reg", "reg": "1000", "part": "all", "deposit": "0",
                      "unit": "3", "output": "deposit", "e-type": "2", "ne": "150",
                      "emin": "0.0", "emax": "1.5"}
    for key, expected in expected_tally.items():
        got = tally.get(key)
        if got is None or got.lower() != expected.lower():
            raise ValueError(f"source input [t-deposit] {key} must be {expected!r}; found {got!r}")


def load_physics_config(path: Path) -> tuple[dict, Path]:
    data = json.loads(path.read_text(encoding="utf-8"))
    expected = {
        "phits_version": "3.341", "region_id": 1000, "region_material": "Pb",
        "dumpall_material_id": 2, "egs5_medium": 2,
        "electron_cutoff_MeV": 0.001, "positron_cutoff_MeV": 0.001,
        "weight_mode": "unit", "energy_unit": "MeV", "deposit_mode": "t-deposit-output-deposit",
    }
    for key, value in expected.items():
        if data.get(key) != value:
            raise ValueError(f"physics config {key} must be {value!r}; found {data.get(key)!r}")
    source_value = data.get("phits_input")
    if not source_value:
        raise ValueError("physics config must name the PHITS input in 'phits_input'")
    source_input = Path(source_value)
    if not source_input.is_absolute():
        source_input = (path.parent / source_input).resolve()
    if not source_input.is_file():
        raise FileNotFoundError(f"PHITS input from physics config not found: {source_input}")
    validate_source_input(source_input)
    return data, source_input


def groups(path: Path):
    stream = path.open("r", encoding="utf-8-sig", newline="")
    reader = csv.DictReader(stream)
    if not reader.fieldnames or "event_id" not in reader.fieldnames:
        stream.close()
        raise ValueError(f"{path} lacks event_id")
    grouped = itertools.groupby(reader, key=lambda row: int(row["event_id"]))
    try:
        for event_id, rows in grouped:
            yield event_id, list(rows)
    finally:
        stream.close()


def take(group_iter, current, event_id):
    if current is None:
        current = next(group_iter, None)
    if current is not None and current[0] < event_id:
        raise ValueError(f"input event {current[0]} is missing from event list")
    if current is not None and current[0] == event_id:
        rows = current[1]
        return rows, next(group_iter, None)
    return [], current


def run(args, writer: EventV2Writer, provenance: dict) -> dict:
    for path in (args.dumpall, args.stopping_table, args.reconstructor, args.track_exporter):
        if not path.is_file():
            raise FileNotFoundError(path)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.work_dir.mkdir(parents=True, exist_ok=True)
    ledger_path = args.work_dir / "physical_ledger.csv"
    event_path = args.work_dir / "raw_events.csv"
    track_path = args.work_dir / "raw_tracks.csv"
    ledger_run = subprocess.run([str(args.reconstructor), str(args.dumpall),
                                 str(args.stopping_table), str(ledger_path)],
                                check=False, text=True, capture_output=True)
    if ledger_run.returncode not in (0, 1) or not ledger_path.is_file():
        raise RuntimeError(f"physical reconstruction failed: rc={ledger_run.returncode}\n{ledger_run.stderr}")
    parse_run = subprocess.run([str(args.track_exporter), str(args.dumpall),
                                 str(event_path), str(track_path)],
                               check=False, text=True, capture_output=True)
    if parse_run.returncode != 0:
        raise RuntimeError(f"raw transport export failed: rc={parse_run.returncode}\n{parse_run.stderr}")

    track_groups = iter(groups(track_path))
    ledger_groups = iter(groups(ledger_path))
    next_tracks = next(track_groups, None)
    next_ledger = next(ledger_groups, None)
    event_count = track_count = particle_count = deposit_count = unresolved_count = 0
    negative_count = duplicate_ids = 0
    with event_path.open("r", encoding="utf-8-sig", newline="") as event_stream:
        for row in csv.DictReader(event_stream):
            event_id = int(row["event_id"])
            if event_id <= 0 or (event_count and event_id <= last_event):
                raise ValueError("raw event IDs must be positive and strictly increasing")
            raw_tracks, next_tracks = take(track_groups, next_tracks, event_id)
            raw_deposits, next_ledger = take(ledger_groups, next_ledger, event_id)

            particles_by_id = {}
            by_source_record = {}
            for track in raw_tracks:
                pid = int(track["particle_id"])
                ityp = int(track["ityp"])
                particle = particles_by_id.setdefault(pid, {"particle_id": pid,
                    "parent_particle_id": None, "ityp": ityp, "tracks": []})
                if particle["ityp"] != ityp:
                    raise ValueError(f"particle {pid} changes ITYP in event {event_id}")
                tid = track["source_record_id"]
                parsed_track = {"source_record_id": tid, "ncol": int(track["ncol"]),
                    "xyz_cm": (float(track["x_cm"]), float(track["y_cm"]), float(track["z_cm"])),
                    "energy_MeV": float(track["energy_MeV"]), "time": float(track["time"]),
                    "region_id": int(track["region_id"])}
                particle["tracks"].append(parsed_track)
                by_source_record[tid] = parsed_track
            particles = list(particles_by_id.values())
            deposits = []
            event_unresolved = 0
            for raw in raw_deposits:
                status = raw["status"]
                energy = float(raw["E_deposit_MeV"]) if raw["E_deposit_MeV"] else None
                negative_count += int(energy is not None and energy < 0)
                kind = "UNRESOLVED" if status != "reconstructed" else "SEGMENT" if raw["position_status"].startswith("segment") else "POINT"
                pos = None
                if kind == "SEGMENT":
                    start = tuple(float(raw[k]) for k in ("x", "y", "z"))
                    end = tuple(float(raw[k]) for k in ("x_end", "y_end", "z_end"))
                    score = None
                    semantics = "EGS5 continuous loss interval from raw Dumpall endpoints; score representative point unavailable"
                elif kind == "POINT" and raw["position_status"] == "secondary_point":
                    pos = tuple(float(raw[k]) for k in ("x", "y", "z"))
                    start = end = None
                    score = pos
                    semantics = "raw child position; local point validated for dead charged e-/e+ secondary branches"
                elif kind == "POINT" and energy == 0.0:
                    record = by_source_record.get(raw["source_line"])
                    if record is None:
                        raise ValueError(f"zero addition {raw['addition_id']} has no raw transport point")
                    pos = record["xyz_cm"]
                    start = end = None
                    score = pos
                    semantics = "raw transport point retained for a zero-energy tally addition; no heat contribution"
                else:
                    start = end = score = None
                    semantics = "unresolved position"
                deposits.append({
                    "addition_id": raw["addition_id"], "particle_id": int(raw["particle_id"]),
                    "process": raw["process"], "region_id": int(raw["region_id"]),
                    "energy_MeV": energy, "kind": kind, "start_xyz_cm": start,
                    "end_xyz_cm": end, "score_xyz_cm": score,
                    "source_record_id": raw["source_line"], "position_semantics": semantics,
                })
                deposit_count += 1
                unresolved_count += int(kind == "UNRESOLVED")
                event_unresolved += int(kind == "UNRESOLVED")
            track_count += len(raw_tracks)
            particle_count += len(particles)
            complete = event_unresolved == 0
            writer.append_event(event_id, particles, deposits, complete=complete,
                                all_absorption="undetermined")
            event_count += 1
            last_event = event_id
    if next_tracks is not None or next_ledger is not None:
        raise ValueError("track or deposit input contains events absent from the raw event stream")
    if negative_count:
        raise ValueError(f"negative physical energies found: {negative_count}")
    if unresolved_count:
        raise ValueError(f"unsupported scoring branches remain ({unresolved_count}); refusing to publish complete HDF5 v2")
    return {
        "input_dumpall": str(args.dumpall), "input_bytes": args.dumpall.stat().st_size,
        "dumpall_sha256": sha256(args.dumpall), "stopping_table_sha256": sha256(args.stopping_table),
        **provenance,
        "physical_ledger_exit_code": ledger_run.returncode,
        "physical_ledger_stdout": ledger_run.stdout.strip(), "transport_export_stdout": parse_run.stdout.strip(),
        "events": event_count, "particles": particle_count, "tracks": track_count,
        "deposits": deposit_count, "unresolved": unresolved_count, "negative": negative_count,
        "output": str(args.output),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dumpall", type=Path, required=True)
    parser.add_argument("--stopping-table", type=Path, required=True)
    parser.add_argument("--source-input", type=Path,
                        help="legacy: PHITS input whose settings are validated")
    parser.add_argument("--physics-config", type=Path,
                        help="explicit JSON contract for the experimental Pb EGS5 mode")
    parser.add_argument("--reconstructor", type=Path, required=True)
    parser.add_argument("--track-exporter", type=Path, required=True)
    parser.add_argument("--converter-executable", type=Path,
                        help="the dump2event executable that launched this driver")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--work-dir", type=Path,
                        help="optional scratch directory (otherwise a private temporary directory is used)")
    args = parser.parse_args()
    try:
        if args.output.exists():
            raise FileExistsError(f"refusing to overwrite existing output: {args.output}")
        if bool(args.source_input) == bool(args.physics_config):
            raise ValueError("provide exactly one of --physics-config or legacy --source-input")
        if args.physics_config:
            physics_data, args.source_input = load_physics_config(args.physics_config)
        else:
            if not args.source_input.is_file():
                raise FileNotFoundError(args.source_input)
            validate_source_input(args.source_input)
            physics_data = None
        # The writer is held open only for this conversion. Output publication is
        # atomic; input data and existing destinations are never modified.
        physics = physics_config(args.source_input)
        provenance = provenance_metadata(physics_data, args.physics_config,
                                         args.source_input, args.dumpall, args.stopping_table,
                                         args.converter_executable, args.reconstructor, args.track_exporter)
        config = {**provenance, "material": "Pb", "medium": "2", "cutoff_MeV": "0.001",
                  "algorithm_version": "pb-egs5-exp-v2-positron",
                  "position_semantics": physics["position_semantics"],
                  "unsupported_processes": physics["unsupported_processes"],
                  "weights": physics["weights"], "diagnostic_mode": True,
                  "tally_semantics": "PHITS [t-deposit] unit=3 output=deposit; transport-step position remains a segment",
                  "complete_semantics": "per-event candidate scoring complete only when no UNRESOLVED rows remain",
                  "dumpall_structure_semantics": "conversion requires terminal NCOL=2",
                  "transport_physics_completeness": "undetermined",
                  "energy_balance_status": "not assessed",
                  "all_absorption": "undetermined",
                  "region_id": physics["region_id"], "dumpall_material_id": physics["dumpall_material_id"],
                  "region_material": physics["region_material"], "egs5_medium": physics["egs5_medium"]}
        started = time.perf_counter()
        with tempfile.TemporaryDirectory(prefix="dump2event-pb-", dir=args.output.parent) as scratch:
            if args.work_dir:
                args.work_dir.mkdir(parents=True, exist_ok=True)
            else:
                args.work_dir = Path(scratch)
            with EventV2Writer(args.output, config) as writer:
                result = run(args, writer, provenance)
        result["output_bytes"] = args.output.stat().st_size
        result["elapsed_seconds"] = time.perf_counter() - started
        print(json.dumps(result, indent=2))
        return 0
    except Exception as error:
        print(f"INVALID: {error}")
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
