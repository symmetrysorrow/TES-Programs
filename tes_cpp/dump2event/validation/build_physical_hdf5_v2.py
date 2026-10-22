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
import os
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


def resolved_esmin_MeV(source_input: Path, phits_version: str) -> float:
    """Resolve PHITS's global [Parameters] esmin, using only the pinned 3.341 default."""
    if phits_version != "3.341":
        raise ValueError("dead-secondary score gate is validated only for PHITS 3.341")
    active = ""
    params: dict[str, str] = {}
    for raw in source_input.read_text(encoding="utf-8-sig", errors="replace").splitlines():
        line = raw.split("#", 1)[0].strip()
        if not line:
            continue
        if line.startswith("[") and line.endswith("]"):
            active = re.sub(r"\s+", "", line[1:-1]).lower()
        elif active == "parameters" and "=" in line:
            key, value = line.split("=", 1)
            params[re.sub(r"\s+", "", key).lower()] = value.strip()
    explicit = params.get("esmin")
    if explicit is None:
        # PHITS 3.341 read02.f maps the `esmin` input name to parz(162).
        # The checked-in default table initializes parz(162) to exactly 1.0e-3 MeV (confirmed from read02.f data parz(161:170), element 2).
        return 1.0e-3
    try:
        value = float(explicit.split()[0])
    except (ValueError, IndexError) as error:
        raise ValueError(f"invalid [Parameters] esmin: {explicit!r}") from error
    if not (value > 0.0 and value < float("inf")):
        raise ValueError("[Parameters] esmin must be finite and positive")
    return value


def material_settings(contract: dict) -> dict:
    """Resolve explicit Pb/Bi single-material contracts without inventing Bi properties."""
    element = contract.get("region_material")
    if element not in ("Pb", "Bi"):
        raise ValueError("region_material must be Pb or Bi (no alloys/other materials)")
    if element == "Pb":
        expected = {"region_id": 1000, "dumpall_material_id": 2, "egs5_medium": 2}
        for key, value in expected.items():
            if contract.get(key) != value:
                raise ValueError(f"validated Pb {key} must be {value}")
        density = contract.get("cell_density_g_cm3", 9.78)
        if type(density) not in (int, float) or float(density) != 9.78:
            raise ValueError("validated Pb cell density must be 9.78 g/cm3")
        return {**expected, "region_material": "Pb", "cell_density_g_cm3": 9.78,
                "validation_status": "PB_SAME_RUN_VALIDATED"}
    for key in ("region_id", "dumpall_material_id", "egs5_medium"):
        value = contract.get(key)
        if type(value) is not int or value <= 0:
            raise ValueError(f"Bi {key} must be an explicitly configured positive integer")
    density = contract.get("cell_density_g_cm3")
    if type(density) not in (int, float) or not (float("inf") > float(density) > 0.0):
        raise ValueError("Bi cell_density_g_cm3 is required; never infer a Pb or tabulated density")
    hash_value = contract.get("stopping_table_sha256")
    if not isinstance(hash_value, str) or not re.fullmatch(r"[0-9a-fA-F]{64}", hash_value):
        raise ValueError("Bi requires the exact private HATCH stopping_table_sha256")
    return {
        "region_material": "Bi", "region_id": contract["region_id"],
        "dumpall_material_id": contract["dumpall_material_id"],
        "egs5_medium": contract["egs5_medium"],
        "cell_density_g_cm3": float(density),
        "stopping_table_sha256": hash_value.lower(),
        "validation_status": "BI_PENDING_PHITS_ORACLE",
    }


def validate_material_stopping_table(contract: dict, table_path: Path) -> str:
    """Verify bytes against the explicitly provided HATCH hash, if available.

    This only establishes file identity and material provenance; it is not
    a PHITS scoring-oracle validation of Bi.
    """
    model = material_settings(contract)
    if not table_path.is_file():
        raise FileNotFoundError(f"EGS5 HATCH table does not exist: {table_path}")
    digest = sha256(table_path)
    expected = model.get("stopping_table_sha256")
    if expected is not None and digest != expected:
        raise ValueError(f"{model['region_material']} HATCH stopping-table SHA-256 mismatch")
    return digest



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
        "esmin_MeV": float(settings["esmin_MeV"]),
        "dead_secondary_score_threshold_MeV": min(float(settings["esmin_MeV"]), 1.0e-6),
        "dead_secondary_score_rule": "PHITS-3.341-tdepstregEVENT-QCLUSTS7-gt-eletb1-v1",
        "weight_mode": str(settings["weight_mode"]).lower(),
        "energy_unit": str(settings["energy_unit"]),
        "deposit_mode": str(settings["deposit_mode"]),
        "algorithm_version": "pb-egs5-exp-v3-score-gate",
        **({"material_validation_status": "BI_PENDING_PHITS_ORACLE",
            "cell_density_g_cm3": material_settings(settings)["cell_density_g_cm3"]}
           if settings["region_material"] == "Bi" else {}),
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
        "esmin_MeV": resolved_esmin_MeV(source_input, "3.341"),
        "dead_secondary_score_threshold_MeV": min(resolved_esmin_MeV(source_input, "3.341"), 1.0e-6),
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
        "esmin_MeV": 1.0e-3,
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
        "reconstruction_algorithm_version": "pb-egs5-exp-v3-score-gate",
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


def validate_source_input(path: Path, contract: dict | None = None) -> None:
    """Validate a one-material PHITS scoring input without asserting Bi oracle closure."""
    model = material_settings(contract or {"region_material": "Pb", "region_id": 1000,
                                                "dumpall_material_id": 2, "egs5_medium": 2})
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
    material_id = model["dumpall_material_id"]
    element = model["region_material"].lower()
    expected_material = rf"mat\s*\[\s*{material_id}\s*\]\s+{element}\s+1(?:\s|$)"
    if not any(re.match(expected_material, line) for line in material_lines):
        raise ValueError(f"source input must define mat[{material_id}] as pure {model['region_material']}")
    cell_lines = sections.get("cell", [])
    region = model["region_id"]
    density = model["cell_density_g_cm3"]
    # Match only the region's single-cell definition, never an unrelated material.
    cell_matches = [
        re.match(rf"^{region}\s+(\d+)\s+(-?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?)(?:\s|$)", line)
        for line in cell_lines if re.match(rf"^{region}(?:\s|$)", line)
    ]
    if len(cell_matches) != 1 or cell_matches[0] is None:
        raise ValueError(f"source input must define exactly one simple cell for region {region}")
    found_material = int(cell_matches[0].group(1))
    found_density = float(cell_matches[0].group(2))
    if found_material != material_id or found_density != -density:
        raise ValueError(f"region {region} requires material {material_id} and density -{density:g}")
    source_definition(sections.get("source", []))
    esmin = resolved_esmin_MeV(path, "3.341")
    # Force evaluation now; the exact score gate is captured in physics config and provenance.
    if not (esmin > 0.0 and esmin < float("inf")):
        raise ValueError("PHITS 3.341 esmin did not resolve to a finite positive value")
    tally = values("t-deposit")
    expected_tally = {"mesh": "reg", "reg": str(model["region_id"]), "part": "all", "deposit": "0",
                      "unit": "3", "output": "deposit", "e-type": "2", "ne": "150",
                      "emin": "0.0", "emax": "1.5"}
    for key, expected in expected_tally.items():
        got = tally.get(key)
        if got is None or got.lower() != expected.lower():
            raise ValueError(f"source input [t-deposit] {key} must be {expected!r}; found {got!r}")


def load_physics_config(path: Path, *, allow_pending_bi: bool = False) -> tuple[dict, Path]:
    data = json.loads(path.read_text(encoding="utf-8"))
    model = material_settings(data)
    expected = {
        "phits_version": "3.341",
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
    validate_source_input(source_input, data)
    if model["region_material"] == "Bi" and not allow_pending_bi:
        raise ValueError("Bi physics is not yet PHITS-oracle validated; HDF5 publication is disabled")
    resolved_esmin = resolved_esmin_MeV(source_input, str(data["phits_version"]))
    if float(data.get("esmin_MeV", float("nan"))) != resolved_esmin:
        raise ValueError(f"physics config esmin_MeV {data.get('esmin_MeV')!r} does not match resolved PHITS value {resolved_esmin}")
    expected_threshold = min(resolved_esmin, 1.0e-6)
    if float(data.get("dead_secondary_score_threshold_MeV", float("nan"))) != expected_threshold:
        raise ValueError("physics config dead_secondary_score_threshold_MeV does not equal min(esmin,1e-6)")
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
    not_scored_path = args.work_dir / "not_scored_candidates.csv"
    event_path = args.work_dir / "raw_events.csv"
    track_path = args.work_dir / "raw_tracks.csv"
    threshold = float(args.dead_secondary_score_threshold_MeV)
    ledger_run = subprocess.run([str(args.reconstructor), str(args.dumpall),
                                 str(args.stopping_table), str(ledger_path),
                                 str(not_scored_path), format(threshold, ".17g")],
                                check=False, text=True, capture_output=True)
    if ledger_run.returncode not in (0, 1) or not ledger_path.is_file():
        raise RuntimeError(f"physical reconstruction failed: rc={ledger_run.returncode}\n{ledger_run.stderr}")
    parse_run = subprocess.run([str(args.track_exporter), str(args.dumpall),
                                 str(event_path), str(track_path)],
                               check=False, text=True, capture_output=True)
    if parse_run.returncode != 0:
        raise RuntimeError(f"raw transport export failed: rc={parse_run.returncode}\n{parse_run.stderr}")
    with not_scored_path.open("r", encoding="utf-8-sig", newline="") as handle:
        not_scored_rows = list(csv.DictReader(handle))
    if any(row.get("classification") != "NOT_SCORED_BY_PHITS" for row in not_scored_rows):
        raise ValueError("non-scored candidate audit contains an unexpected classification")

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
        "not_scored_by_phits_candidates": len(not_scored_rows),
        "not_scored_candidate_audit": str(not_scored_path),
        "events": event_count, "particles": particle_count, "tracks": track_count,
        "deposits": deposit_count, "unresolved": unresolved_count, "negative": negative_count,
        "output": str(args.output),
    }


# Pinned to the exact, same-run-PHITS-validated historical HATCH CSV.
# This data file is NOT automatically published here: verify distribution
# permission before adding it to a public Git repository.
BUNDLED_PB_TABLE_NAME = "pb_egs5_phits3341.csv"
BUNDLED_PB_TABLE_SHA256 = "20fd86eb0b57483bd19f20a8388d2c1ecb0431a02fa99f0933e2c3c4748514e9"


def resolve_stopping_table(explicit: Path | None) -> Path:
    """Select an explicitly provided table or an exact, checksum-pinned bundled Pb table.

    No unverified physics table is silently substituted, and no synthetic
    fixture is ever a valid default.
    """
    if explicit is not None:
        candidates = [explicit]
    else:
        # Prefer the immutable, SHA-pinned distribution resource over any old
        # per-user environment variable. The environment path is compatibility
        # fallback only when no bundled table is installed.
        home = Path(__file__).resolve().parent
        candidates = [
            home / BUNDLED_PB_TABLE_NAME,  # Python package install beside driver
            home.parent / "data" / BUNDLED_PB_TABLE_NAME,  # repository source tree
            home.parent.parent / "bin" / BUNDLED_PB_TABLE_NAME,  # CMake prefix/bin
        ]
        configured = os.environ.get("TES_CPP_PB_HATCH_TABLE")
        if configured:
            candidates.append(Path(configured).expanduser())
    for candidate in candidates:
        if candidate.is_file():
            actual = sha256(candidate)
            if actual != BUNDLED_PB_TABLE_SHA256:
                raise ValueError(
                    f"Pb HATCH coefficient checksum mismatch: {candidate} "
                    f"(expected {BUNDLED_PB_TABLE_SHA256}, got {actual})"
                )
            return candidate
    if explicit is not None:
        raise FileNotFoundError(f"explicit stopping table was not found: {explicit}")
    raise FileNotFoundError(
        "verified Pb HATCH coefficient table is not bundled. "
        "Install the byte-verified table as tes_cpp/dump2event/data/"
        + BUNDLED_PB_TABLE_NAME +
        "; set TES_CPP_PB_HATCH_TABLE for a private installation, or specify "
        "--stopping-table with the exact SHA-256-pinned CSV."
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dumpall", type=Path, required=True)
    parser.add_argument("--stopping-table", type=Path,
                        help="optional override; omitted means the checksum-pinned bundled Pb HATCH CSV")
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
        args.stopping_table = resolve_stopping_table(args.stopping_table)
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
        args.dead_secondary_score_threshold_MeV = float(
            physics_data["dead_secondary_score_threshold_MeV"] if physics_data else
            physics["dead_secondary_score_threshold_MeV"])
        provenance = provenance_metadata(physics_data, args.physics_config,
                                         args.source_input, args.dumpall, args.stopping_table,
                                         args.converter_executable, args.reconstructor, args.track_exporter)
        config = {**provenance, "material": physics_data["region_material"],
                  "medium": str(physics_data["egs5_medium"]), "cutoff_MeV": "0.001",
                  "algorithm_version": "pb-egs5-exp-v3-score-gate",
                  "position_semantics": physics["position_semantics"],
                  "unsupported_processes": physics["unsupported_processes"],
                  "weights": physics["weights"], "diagnostic_mode": True,
                  "tally_semantics": "PHITS [t-deposit] unit=3 output=deposit; transport-step position remains a segment",
                  "complete_semantics": "per-event scoring complete only when no in-scope scored candidate is UNRESOLVED; proven not-scored candidates are excluded from deposits and preserved in audit",
                  "dead_secondary_score_rule": "PHITS-3.341-tdepstregEVENT-QCLUSTS7-gt-eletb1-v1",
                  "dead_secondary_score_threshold_MeV": format(float(physics_data["dead_secondary_score_threshold_MeV"] if physics_data else physics["dead_secondary_score_threshold_MeV"]), ".17g"),
                  "esmin_MeV": format(float(physics_data["esmin_MeV"] if physics_data else physics["esmin_MeV"]), ".17g"),
                  "dumpall_structure_semantics": "conversion requires terminal NCOL=2",
                  "transport_physics_completeness": "undetermined",
                  "energy_balance_status": "not assessed",
                  "all_absorption": "undetermined",
                  "region_id": physics["region_id"], "dumpall_material_id": physics["dumpall_material_id"],
                  "region_material": physics["region_material"], "egs5_medium": physics["egs5_medium"]}
        started = time.perf_counter()
        # load_physics_config refuses Bi until an independent same-run PHITS
        # scoring oracle validates that material and source configuration.
        with tempfile.TemporaryDirectory(prefix="dump2event-egs5-", dir=args.output.parent) as scratch:
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
