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
import subprocess
import time
from pathlib import Path

from tes_cpp.event_hdf5_v2 import EventV2Writer, iter_deposits, iter_events


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def physics_config(source_input: Path):
    return {
        "phits_version": "3.341", "region_id": 1000, "region_material": "PB",
        "dumpall_material_id": 2, "egs5_medium": 2, "electron_cutoff_MeV": 0.001,
        "weights": "unit-weight validated slice only", "diagnostic_mode": True,
        "unsupported_processes": "stored UNRESOLVED; never old-formula fallback",
        "source_input_sha256": sha256(source_input),
        "position_semantics": "raw transport chords and raw child point; scorer representative coordinate kept separate",
    }


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


def run(args, writer: EventV2Writer) -> dict:
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

    config = physics_config(args.source_input)
    config_hash = hashlib.sha256(json.dumps(config, sort_keys=True, separators=(",", ":")).encode()).hexdigest()
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
                    semantics = "raw child position; local point validated only for the supported dead-electron branch"
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
    return {
        "input_dumpall": str(args.dumpall), "input_bytes": args.dumpall.stat().st_size,
        "dumpall_sha256": sha256(args.dumpall), "stopping_table_sha256": sha256(args.stopping_table),
        "source_input_sha256": config["source_input_sha256"],
        "physics_config_sha256": config_hash, "phits_version": "3.341",
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
    parser.add_argument("--source-input", type=Path, required=True,
                        help="PHITS input whose settings are summarized in the provenance")
    parser.add_argument("--reconstructor", type=Path, required=True)
    parser.add_argument("--track-exporter", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--work-dir", type=Path, required=True)
    args = parser.parse_args()
    try:
        if args.output.exists():
            raise FileExistsError(f"refusing to overwrite existing output: {args.output}")
        if not args.source_input.is_file():
            raise FileNotFoundError(args.source_input)
        # The writer is held open only for this conversion. Output publication is
        # atomic; input data and existing destinations are never modified.
        physics = physics_config(args.source_input)
        config = {"phits_version": "3.341", "dumpall_sha256": sha256(args.dumpall),
                  "physics_config_sha256": hashlib.sha256(json.dumps(physics, sort_keys=True, separators=(",", ":")).encode()).hexdigest(),
                  "stopping_table_sha256": sha256(args.stopping_table), "material": "Pb",
                  "medium": "2", "cutoff_MeV": "0.001", "algorithm_version": "pb-egs5-exp-v1",
                  "position_semantics": physics["position_semantics"],
                  "source_input_sha256": physics["source_input_sha256"],
                  "unsupported_processes": physics["unsupported_processes"],
                  "weights": physics["weights"], "diagnostic_mode": True,
                  "region_id": physics["region_id"], "dumpall_material_id": physics["dumpall_material_id"],
                  "region_material": physics["region_material"], "egs5_medium": physics["egs5_medium"]}
        started = time.perf_counter()
        with EventV2Writer(args.output, config) as writer:
            result = run(args, writer)
        result["output_bytes"] = args.output.stat().st_size
        result["elapsed_seconds"] = time.perf_counter() - started
        print(json.dumps(result, indent=2))
        return 0
    except Exception as error:
        print(f"INVALID: {error}")
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
