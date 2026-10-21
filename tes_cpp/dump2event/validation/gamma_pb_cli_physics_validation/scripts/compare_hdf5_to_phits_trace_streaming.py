#!/usr/bin/env python3
"""Bounded-memory same-run HDF5/Dumpall/PHITS trace comparison."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path

import h5py


def sha(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def rounded(value: str | float) -> float:
    return round(float(value), 13)


def raw_key(row: dict) -> tuple:
    return (int(row["event_id"]), int(row["particle_id"]), int(row["ityp"]),
            int(row["ncol"]), rounded(row["energy_MeV"]), rounded(row["x_cm"]),
            rounded(row["y_cm"]), rounded(row["z_cm"]), int(row["region_id"]))


def raw_sequence_map(run_dir: Path, tracks_csv: Path):
    """Join ordered transport rows to ordered DUMP sidecar records without a map."""
    tracks_stream = tracks_csv.open(newline="", encoding="utf-8-sig")
    tracks = iter(csv.DictReader(tracks_stream))
    current = next(tracks, None)
    try:
        with (run_dir / "dump_record_sequence.log").open(encoding="ascii") as sidecar:
            for line_no, line in enumerate(sidecar, 1):
                fields = line.split()
                if not fields or fields[0] != "DUMP":
                    continue
                if len(fields) < 23:
                    raise ValueError(f"short DUMP sidecar row at line {line_no}")
                ncol = int(fields[7])
                if not (10 <= ncol <= 15):
                    continue
                if current is None:
                    raise ValueError("DUMP sidecar has a transport record missing from raw_tracks.csv")
                side = {"event_id": fields[2], "particle_id": fields[5], "ityp": fields[6],
                        "ncol": fields[7], "energy_MeV": fields[13], "x_cm": fields[15],
                        "y_cm": fields[16], "z_cm": fields[17], "region_id": fields[21]}
                if raw_key(side) != raw_key(current):
                    raise ValueError(f"ordered raw/DUMP sidecar identity mismatch at DUMP sequence {fields[1]}")
                yield int(current["source_record_id"]), int(fields[1]), current
                current = next(tracks, None)
        if current is not None:
            raise ValueError("raw_tracks.csv contains transport rows absent from DUMP sequence sidecar")
    finally:
        tracks_stream.close()


def internal_rows(run_dir: Path, internal_csv: Path):
    """Pair TALLYSEQ and PRE/POST rows in stream order and expose each increment."""
    with internal_csv.open(newline="", encoding="utf-8-sig") as stream:
        pres = iter(csv.DictReader(stream))
        previous = -1
        with (run_dir / "tally_dumpseq.log").open(encoding="ascii") as sidecar:
            for line_no, line in enumerate(sidecar, 1):
                fields = line.split()
                if not fields or fields[0] != "TALLYSEQ":
                    continue
                if len(fields) < 17:
                    raise ValueError(f"short TALLYSEQ sidecar row at line {line_no}")
                pre = next(pres, None)
                if pre is None:
                    raise ValueError("PHITS TALLYSEQ has more rows than paired PRE/POST trace")
                serial = int(fields[1])
                if serial < previous:
                    raise ValueError("TALLYSEQ serials are not nondecreasing")
                previous = serial
                identity = {"event_id": int(float(fields[2])), "particle_id": int(fields[5]),
                            "ityp": int(fields[6]), "jcoll": int(fields[7]), "ncol": int(fields[8]),
                            "region_id": int(fields[10]), "weight": float(fields[13]),
                            "E_MeV": float(fields[14]), "EC_MeV": float(fields[15]),
                            "tlv_MeV": float(fields[16])}
                for key, value in (("nocas", identity["event_id"]), ("particle_id", identity["particle_id"]),
                                   ("ityp", identity["ityp"]), ("jcoll", identity["jcoll"]),
                                   ("ncol", identity["ncol"])):
                    if int(float(pre[key])) != value:
                        raise ValueError(f"PRE/TALLYSEQ {key} mismatch at TALLYSEQ line {line_no}")
                increment = float(pre["tally_increment_MeV"])
                if abs(increment - identity["tlv_MeV"]) > 1e-12:
                    raise ValueError(f"PRE/POST increment differs from tlv at TALLYSEQ line {line_no}")
                yield serial, identity, pre
        if next(pres, None) is not None:
            raise ValueError("paired PRE/POST trace has more rows than PHITS TALLYSEQ")


def text_value(value):
    return value.decode("utf-8") if isinstance(value, bytes) else str(value)


def vec(value):
    return tuple(float(v) for v in value)


def max_delta(left, right):
    return max(abs(a - b) for a, b in zip(left, right))


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-dir", type=Path, required=True)
    parser.add_argument("--dumpall", type=Path, required=True,
                        help="the exact raw Dumpall file used for this PHITS run")
    parser.add_argument("--tracks-csv", type=Path, required=True)
    parser.add_argument("--internal-csv", type=Path, required=True)
    parser.add_argument("--hdf5", type=Path, required=True)
    parser.add_argument("--output-json", type=Path, required=True)
    parser.add_argument("--output-csv", type=Path)
    args = parser.parse_args()

    raw_rows = iter(raw_sequence_map(args.run_dir, args.tracks_csv))
    raw = next(raw_rows, None)
    current_source_line = None
    internals = iter(internal_rows(args.run_dir, args.internal_csv))
    internal = next(internals, None)
    counts = Counter()
    energy_max = energy_sum_abs = increment_max = event_max = 0.0
    energy_diffs = []
    position_max = 0.0
    event_h5 = defaultdict(float)
    event_tlv = defaultdict(float)
    event_counts = Counter()
    output = None
    writer = None
    if args.output_csv:
        args.output_csv.parent.mkdir(parents=True, exist_ok=True)
        output = args.output_csv.open("w", newline="", encoding="utf-8")
        writer = csv.writer(output)
        writer.writerow(("addition_id", "source_record_id", "dump_sequence", "event_id", "particle_id",
                         "region_hdf5", "region_internal", "energy_hdf5_MeV", "energy_tlv_MeV",
                         "abs_diff_MeV", "weight", "position_kind", "position_max_abs_diff_cm"))

    try:
        with h5py.File(args.hdf5, "r") as h5:
            d = h5["deposits"]
            total = len(d["event_id"])
            selected = ("addition_id", "source_record_id", "event_id", "particle_id", "region_id",
                        "energy_MeV", "energy_valid", "kind", "start_position_valid",
                        "end_position_valid", "score_position_valid", "start_xyz_cm", "end_xyz_cm", "score_xyz_cm")
            for base in range(0, total, 32768):
                stop = min(base + 32768, total)
                block = {name: d[name][base:stop] for name in selected}
                for index in range(stop - base):
                    counts["hdf5_rows"] += 1
                    sid = int(text_value(block["source_record_id"][index]))
                    if sid != current_source_line:
                        if current_source_line is not None and raw is not None and raw[0] == current_source_line:
                            raw = next(raw_rows, None)
                        while raw is not None and raw[0] < sid:
                            counts["raw_rows_not_scored"] += 1
                            raw = next(raw_rows, None)
                        current_source_line = sid
                    if raw is None or raw[0] != sid:
                        counts["unmapped_hdf5_source_line"] += 1
                        raise ValueError(f"HDF5 source line {sid} has no ordered raw-record sequence")
                    serial, raw_track = raw[1], raw[2]
                    if internal is None:
                        counts["missing_internal"] += 1
                        raise ValueError("PHITS internal tally trace ended before HDF5 additions")
                    while internal[0] < serial:
                        counts["internal_not_reconstructed"] += 1
                        internal = next(internals, None)
                        if internal is None:
                            raise ValueError("PHITS internal tally rows ended before HDF5 additions")
                    if internal[0] != serial:
                        counts["missing_internal"] += 1
                        raise ValueError(f"Dump sequence {serial} has no matching PHITS internal tally addition")
                    identity = internal[1]
                    pre = internal[2]
                    event = int(block["event_id"][index])
                    particle = int(block["particle_id"][index])
                    region = int(block["region_id"][index])
                    if event != identity["event_id"]: counts["event_mismatch"] += 1
                    if particle != identity["particle_id"]: counts["particle_mismatch"] += 1
                    if region != identity["region_id"]: counts["region_mismatch"] += 1
                    if not bool(block["energy_valid"][index]):
                        counts["invalid_energy"] += 1
                        raise ValueError("HDF5 addition has invalid energy")
                    h_energy = float(block["energy_MeV"][index])
                    delta = abs(h_energy - identity["tlv_MeV"])
                    energy_max = max(energy_max, delta)
                    energy_sum_abs += delta
                    energy_diffs.append(delta)
                    event_h5[event] += h_energy
                    event_tlv[event] += identity["tlv_MeV"]
                    event_counts[event] += 1
                    kind = int(block["kind"][index])
                    position_delta = None
                    p0 = tuple(float(pre[k]) for k in ("x_start_cm", "y_start_cm", "z_start_cm"))
                    p1 = tuple(float(pre[k]) for k in ("x_end_cm", "y_end_cm", "z_end_cm"))
                    if kind == 2 and bool(block["start_position_valid"][index]) and bool(block["end_position_valid"][index]):
                        position_delta = max(max_delta(vec(block["start_xyz_cm"][index]), p0),
                                             max_delta(vec(block["end_xyz_cm"][index]), p1))
                    elif bool(block["score_position_valid"][index]):
                        point = vec(block["score_xyz_cm"][index])
                        position_delta = min(max_delta(point, p0), max_delta(point, p1))
                    if position_delta is None:
                        counts["position_unavailable"] += 1
                    else:
                        counts["position_rows_compared"] += 1
                        position_max = max(position_max, position_delta)
                    if writer:
                        writer.writerow((text_value(block["addition_id"][index]), sid, serial, event, particle,
                                         region, identity["region_id"], format(h_energy, ".17g"),
                                         format(identity["tlv_MeV"], ".17g"), format(delta, ".17g"),
                                         format(identity["weight"], ".17g"), kind,
                                         "" if position_delta is None else format(position_delta, ".17g")))
                    increment_max = max(increment_max, abs(float(pre["tally_increment_MeV"]) - identity["tlv_MeV"]))
                    counts["matched_rows"] += 1
                    internal = next(internals, None)

        event_keys = set(event_h5) | set(event_tlv)
        if event_keys:
            event_max = max(abs(event_h5[event] - event_tlv[event]) for event in event_keys)
        while internal is not None:
            counts["internal_not_reconstructed"] += 1
            internal = next(internals, None)
        if raw is not None and current_source_line is not None and raw[0] == current_source_line:
            raw = next(raw_rows, None)
        while raw is not None:
            counts["raw_rows_not_scored"] += 1
            raw = next(raw_rows, None)
        if counts["matched_rows"] != counts["hdf5_rows"]:
            raise ValueError("not every HDF5 addition matched an internal tally row")
        result = {
            **counts,
            "internal_rows_not_reconstructed": counts["internal_not_reconstructed"],
            "missing_internal": counts["missing_internal"],
            "event_mismatches": counts["event_mismatch"],
            "particle_mismatches": counts["particle_mismatch"],
            "region_mismatches": counts["region_mismatch"],
            "position_unavailable": counts["position_unavailable"],
            "event_count": len(event_keys),
            "hdf5_vs_tlv_max_abs_MeV": energy_max,
            "hdf5_vs_tlv_sum_abs_MeV": energy_sum_abs,
            "hdf5_vs_tlv_rows_over_1e-12_MeV": sum(delta > 1e-12 for delta in energy_diffs),
            "pre_increment_vs_tlv_max_abs_MeV": increment_max,
            "hdf5_total_MeV": sum(event_h5.values()),
            "internal_tlv_total_MeV": sum(event_tlv.values()),
            "event_sum_max_abs_diff_MeV": event_max,
            "position_max_abs_diff_cm": position_max,
        "dumpall_sha256": sha(args.dumpall),
        "dumpall_bytes": args.dumpall.stat().st_size,
            "hdf5_sha256": sha(args.hdf5),
            "internal_trace_sha256": sha(args.internal_csv),
            "tally_dumpseq_sha256": sha(args.run_dir / "tally_dumpseq.log"),
        }
        args.output_json.parent.mkdir(parents=True, exist_ok=True)
        args.output_json.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
        print(json.dumps(result, indent=2))
    finally:
        if output:
            output.close()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
