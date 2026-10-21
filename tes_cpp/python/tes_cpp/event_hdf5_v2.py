"""Appendable flat HDF5 v2 schema for PHITS transport and deposit records.

This module deliberately keeps transport segments distinct from point scores.
It is a storage/API layer; it does not reconstruct physics or assign TES blocks.
"""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import re
import tempfile
from typing import Mapping, Sequence

import h5py
import numpy as np

FORMAT_NAME = "tes-dump2event"
FORMAT_VERSION = 2
DEPOSIT_KINDS = {"POINT": 1, "SEGMENT": 2, "UNRESOLVED": 3}
ABSORPTION_STATES = {"confirmed", "rejected", "undetermined"}
REQUIRED_PROVENANCE = (
    "phits_version", "dumpall_sha256", "physics_config_sha256",
    "stopping_table_sha256", "material", "medium", "cutoff_MeV",
    "algorithm_version", "position_semantics",
)


def _text(value) -> str:
    if isinstance(value, bytes):
        return value.decode("utf-8")
    return str(value)


def _append(dataset, values) -> tuple[int, int]:
    start = len(dataset)
    array = np.asarray(values, dtype=dataset.dtype)
    if array.ndim == 0:
        array = array.reshape(1)
    dataset.resize(start + len(array), axis=0)
    if len(array):
        dataset[start:] = array
    return start, start + len(array)


def _vector(value, field: str, *, required: bool = False):
    if value is None:
        if required:
            raise ValueError(f"deposit is missing {field}")
        return (0.0, 0.0, 0.0)
    if len(value) != 3:
        raise ValueError(f"{field} must contain exactly 3 coordinates")
    result = tuple(float(component) for component in value)
    if not np.isfinite(result).all():
        raise ValueError(f"{field} must contain finite coordinates")
    return result


class EventV2Writer:
    """Stream events to a temporary file, publishing only after validation.

    ``append_event`` accepts flat particle/track/deposit mappings. An event with
    unresolved deposits must be marked incomplete; unknown energy is represented
    by ``energy_valid=False`` and is never silently interpreted as a known zero.
    """

    def __init__(self, path, provenance: Mapping[str, object], *, compression="gzip"):
        self.path = Path(path)
        missing = [name for name in REQUIRED_PROVENANCE if name not in provenance or provenance[name] in (None, "")]
        if missing:
            raise ValueError(f"missing required provenance: {', '.join(missing)}")
        self.path.parent.mkdir(parents=True, exist_ok=True)
        fd, temporary = tempfile.mkstemp(prefix=f".{self.path.name}.", suffix=".tmp", dir=self.path.parent)
        os.close(fd)
        self.temporary_path = Path(temporary)
        self.file = h5py.File(self.temporary_path, "w")
        self.file.attrs["format"] = FORMAT_NAME
        self.file.attrs["format_version"] = FORMAT_VERSION
        self.file.attrs["energy_unit"] = "MeV"
        self.file.attrs["length_unit"] = "cm"
        self.file.attrs["schema"] = "flat-offset-v2"
        self.provenance = self.file.create_group("provenance")
        for key, value in provenance.items():
            if key == "format":
                continue
            if key in REQUIRED_PROVENANCE:
                self.provenance.attrs[key] = str(value)
            else:
                self.provenance.attrs[key] = value if isinstance(
                    value, (str, int, float, bool, np.integer, np.floating)) else str(value)
        hash_keys = ["dumpall_sha256", "physics_config_sha256", "stopping_table_sha256"]
        for key in ("physics_contract_sha256", "resolved_physics_sha256", "phits_input_sha256",
                    "converter_executable_sha256", "reconstructor_executable_sha256",
                    "transport_exporter_sha256", "driver_sha256", "phits_executable_sha256"):
            if key in provenance:
                hash_keys.append(key)
        for key in hash_keys:
            if not re.fullmatch(r"[0-9a-fA-F]{64}", str(provenance[key])):
                self.file.close()
                self.temporary_path.unlink(missing_ok=True)
                raise ValueError(f"{key} must be a 64-character SHA-256 hex digest")
        revision = int(provenance.get("provenance_schema_revision", 1))
        if revision >= 2:
            if ("resolved_physics_sha256" not in provenance or
                    "physics_contract_sha256" not in provenance or
                    "resolved_physics_json" not in provenance):
                self.file.close()
                self.temporary_path.unlink(missing_ok=True)
                raise ValueError("provenance revision 2 requires contract SHA, resolved JSON, and resolved SHA")
            if str(provenance["physics_config_sha256"]).lower() != str(provenance["resolved_physics_sha256"]).lower():
                self.file.close()
                self.temporary_path.unlink(missing_ok=True)
                raise ValueError("physics_config_sha256 is a compatibility alias and must equal resolved_physics_sha256")
            try:
                parsed = json.loads(str(provenance["resolved_physics_json"]))
                canonical = json.dumps(parsed, sort_keys=True, separators=(",", ":"), ensure_ascii=False).encode("utf-8")
            except (TypeError, ValueError) as error:
                self.file.close()
                self.temporary_path.unlink(missing_ok=True)
                raise ValueError("resolved_physics_json must be valid canonical JSON") from error
            if canonical.decode("utf-8") != str(provenance["resolved_physics_json"]):
                self.file.close()
                self.temporary_path.unlink(missing_ok=True)
                raise ValueError("resolved_physics_json must use sorted compact canonical UTF-8 JSON")
            resolved_digest = hashlib.sha256(canonical).hexdigest()
            if resolved_digest != str(provenance["resolved_physics_sha256"]).lower():
                self.file.close()
                self.temporary_path.unlink(missing_ok=True)
                raise ValueError("resolved_physics_sha256 does not match resolved_physics_json")

        self.events = self.file.create_group("events")
        self.particles = self.file.create_group("particles")
        self.tracks = self.file.create_group("tracks")
        self.deposits = self.file.create_group("deposits")
        self._event_ids = set()
        self._addition_ids = set()
        self._datasets(self.events, {
            "event_id": np.int64, "particle_offsets": np.int64,
            "deposit_offsets": np.int64, "deposit_total_MeV": np.float64,
            "complete": np.bool_, "all_absorption": h5py.string_dtype("utf-8"),
            "unresolved_count": np.uint64,
        }, compression)
        self.events["particle_offsets"].resize(1, axis=0)
        self.events["particle_offsets"][0] = 0
        self.events["deposit_offsets"].resize(1, axis=0)
        self.events["deposit_offsets"][0] = 0
        self._datasets(self.particles, {
            "event_id": np.int64, "particle_id": np.int64,
            "parent_particle_id": np.int64, "parent_valid": np.bool_,
            "ityp": np.int32, "track_offsets": np.int64,
        }, compression)
        self.particles["track_offsets"].resize(1, axis=0)
        self.particles["track_offsets"][0] = 0
        self._datasets(self.tracks, {
            "source_record_id": h5py.string_dtype("utf-8"), "ncol": np.int16,
            "particle_id": np.int64, "x_cm": np.float64, "y_cm": np.float64,
            "z_cm": np.float64, "energy_MeV": np.float64, "time": np.float64,
            "region_id": np.int32, "time_valid": np.bool_,
        }, compression)
        self._datasets(self.deposits, {
            "addition_id": h5py.string_dtype("utf-8"), "event_id": np.int64,
            "particle_id": np.int64, "process": h5py.string_dtype("utf-8"),
            "region_id": np.int32, "energy_MeV": np.float64,
            "energy_valid": np.bool_, "kind": np.uint8,
            "start_xyz_cm": np.float64, "end_xyz_cm": np.float64,
            "start_position_valid": np.bool_, "end_position_valid": np.bool_,
            "score_xyz_cm": np.float64, "score_position_valid": np.bool_,
            "source_record_id": h5py.string_dtype("utf-8"),
            "position_semantics": h5py.string_dtype("utf-8"),
        }, compression)
        self._buffers = {
            group_name: {name: [] for name in group.keys()}
            for group_name, group in (("events", self.events), ("particles", self.particles),
                                      ("tracks", self.tracks), ("deposits", self.deposits))
        }
        self._buffered_events = 0
        self._failed = False

    def _total(self, group: str, name: str) -> int:
        return len(getattr(self, group)[name]) + len(self._buffers[group][name])

    def _push(self, group: str, name: str, values) -> None:
        self._buffers[group][name].extend(values)

    def _flush(self) -> None:
        try:
            for group_name, fields in self._buffers.items():
                group = getattr(self, group_name)
                for name, values in fields.items():
                    if values:
                        _append(group[name], values)
                        fields[name] = []
            self._buffered_events = 0
            self.file.flush()
        except Exception:
            # HDF5 dataset writes cannot be rolled back as one transaction.
            # Poison this private temporary writer and forbid publication or
            # further appends after a storage failure.
            self._failed = True
            raise

    @staticmethod
    def _datasets(group, definitions, compression):
        for name, dtype in definitions.items():
            kwargs = {"shape": (0,), "maxshape": (None,), "dtype": dtype, "chunks": (1024,)}
            if compression:
                kwargs.update(compression=compression, shuffle=True if dtype not in (h5py.string_dtype("utf-8"),) else False)
            if name in ("start_xyz_cm", "end_xyz_cm", "score_xyz_cm"):
                kwargs["shape"] = (0, 3)
                kwargs["maxshape"] = (None, 3)
                kwargs["chunks"] = (1024, 3)
            group.create_dataset(name, **kwargs)

    def append_event(self, event_id: int, particles: Sequence[Mapping], deposits: Sequence[Mapping], *,
                     complete: bool, all_absorption: str = "undetermined"):
        if self.file is None:
            raise RuntimeError("writer is closed")
        if self._failed:
            raise RuntimeError("writer is unusable after a storage failure")
        event_id = int(event_id)
        if event_id < 0 or event_id > np.iinfo(np.int64).max:
            raise ValueError("event_id is outside the HDF5 int64 range")
        if all_absorption not in ABSORPTION_STATES:
            raise ValueError(f"invalid all_absorption state: {all_absorption}")
        if all_absorption == "confirmed" and not complete:
            raise ValueError("all_absorption=confirmed requires a complete event")
        if event_id in self._event_ids:
            raise ValueError(f"duplicate event_id {event_id}")
        pstart = self._total("particles", "particle_id")
        dstart = self._total("deposits", "event_id")
        tstart = self._total("tracks", "ncol")
        track_count = sum(len(particle.get("tracks", ())) for particle in particles)
        for label, value in (("particle offset", pstart + len(particles)),
                             ("deposit offset", dstart + len(deposits)),
                             ("track offset", tstart + track_count)):
            if value > np.iinfo(np.int64).max:
                raise OverflowError(f"{label} exceeds the HDF5 int64 range")
        particle_ids = set()
        new_addition_ids = set()
        known_total = 0.0
        unresolved = 0
        pbuf = {name: [] for name in ("event_id", "particle_id", "parent_particle_id", "parent_valid", "ityp")}
        tbuf = {name: [] for name in ("source_record_id", "ncol", "particle_id", "x_cm", "y_cm", "z_cm",
                                      "energy_MeV", "time", "region_id", "time_valid")}
        dbuf = {name: [] for name in ("addition_id", "event_id", "particle_id", "process", "region_id",
                                      "energy_MeV", "energy_valid", "kind", "start_xyz_cm", "end_xyz_cm",
                                      "start_position_valid", "end_position_valid", "score_xyz_cm",
                                      "score_position_valid", "source_record_id", "position_semantics")}
        track_offsets = []
        for particle in particles:
            pid = int(particle["particle_id"])
            if pid < 0 or pid > np.iinfo(np.int64).max:
                raise ValueError("particle_id is outside the HDF5 int64 range")
            if pid in particle_ids:
                raise ValueError(f"duplicate particle_id {pid} in event {event_id}")
            particle_ids.add(pid)
            has_parent = particle.get("parent_particle_id") is not None
            pbuf["event_id"].append(event_id)
            pbuf["particle_id"].append(pid)
            pbuf["parent_particle_id"].append(int(particle["parent_particle_id"]) if has_parent else -1)
            pbuf["parent_valid"].append(has_parent)
            pbuf["ityp"].append(int(particle["ityp"]))
            for track in particle.get("tracks", ()):
                xyz = _vector(track.get("xyz_cm"), "track.xyz_cm", required=True)
                energy = float(track["energy_MeV"])
                if not np.isfinite(energy) or energy < 0.0:
                    raise ValueError("track energy must be finite and nonnegative")
                source_record_id = str(track["source_record_id"])
                if not source_record_id:
                    raise ValueError("track source_record_id must be nonempty")
                ncol = int(track["ncol"])
                if ncol < 0 or ncol > np.iinfo(np.int16).max:
                    raise ValueError("track ncol is outside the HDF5 int16 range")
                region_id = int(track["region_id"])
                if region_id < 0 or region_id > np.iinfo(np.int32).max:
                    raise ValueError("track region_id is outside the HDF5 int32 range")
                tbuf["source_record_id"].append(source_record_id)
                tbuf["ncol"].append(ncol)
                tbuf["particle_id"].append(pid)
                tbuf["energy_MeV"].append(energy)
                for axis, value in zip(("x_cm", "y_cm", "z_cm"), xyz):
                    tbuf[axis].append(value)
                time_valid = track.get("time") is not None
                time = float(track["time"]) if time_valid else 0.0
                if not np.isfinite(time):
                    raise ValueError("track time must be finite")
                tbuf["time"].append(time)
                tbuf["time_valid"].append(time_valid)
                tbuf["region_id"].append(region_id)
            track_offsets.append(tstart + len(tbuf["ncol"]))

        for particle in particles:
            parent = particle.get("parent_particle_id")
            if parent is not None and int(parent) not in particle_ids:
                raise ValueError(f"parent_particle_id {parent} is absent from event {event_id}")

        for deposit in deposits:
            addition_id = str(deposit["addition_id"])
            if not addition_id:
                raise ValueError("addition_id must be nonempty")
            if addition_id in self._addition_ids or addition_id in new_addition_ids:
                raise ValueError(f"duplicate addition_id {addition_id}")
            new_addition_ids.add(addition_id)
            kind_name = str(deposit["kind"]).upper()
            if kind_name not in DEPOSIT_KINDS:
                raise ValueError(f"invalid deposit kind {kind_name!r}")
            kind = DEPOSIT_KINDS[kind_name]
            valid = deposit.get("energy_MeV") is not None
            energy = float(deposit["energy_MeV"]) if valid else 0.0
            if valid and (not np.isfinite(energy) or energy < 0.0):
                raise ValueError("known physical deposit energy must be finite and nonnegative")
            if kind == DEPOSIT_KINDS["UNRESOLVED"] and valid:
                raise ValueError("UNRESOLVED deposits must use energy_MeV=None and energy_valid=False")
            if kind != DEPOSIT_KINDS["UNRESOLVED"] and not valid:
                raise ValueError(f"{kind_name} deposit requires a known energy")
            if kind == DEPOSIT_KINDS["SEGMENT"]:
                start = _vector(deposit.get("start_xyz_cm"), "start_xyz_cm", required=True)
                end = _vector(deposit.get("end_xyz_cm"), "end_xyz_cm", required=True)
            else:
                start = _vector(deposit.get("start_xyz_cm"), "start_xyz_cm")
                end = _vector(deposit.get("end_xyz_cm"), "end_xyz_cm")
            start_valid = deposit.get("start_xyz_cm") is not None
            end_valid = deposit.get("end_xyz_cm") is not None
            score_valid = deposit.get("score_xyz_cm") is not None
            if kind == DEPOSIT_KINDS["POINT"] and not score_valid:
                raise ValueError("POINT deposit requires score_xyz_cm")
            score = _vector(deposit.get("score_xyz_cm"), "score_xyz_cm")
            pid = int(deposit.get("particle_id", -1))
            if pid >= 0 and pid not in particle_ids:
                raise ValueError(f"deposit references absent particle_id {pid} in event {event_id}")
            if pid < -1:
                raise ValueError("deposit particle_id must be -1 or nonnegative")
            process = str(deposit["process"])
            source_record_id = str(deposit["source_record_id"])
            if not process or not source_record_id:
                raise ValueError("deposit process and source_record_id must be nonempty")
            region_id = int(deposit["region_id"])
            if region_id < 0 or region_id > np.iinfo(np.int32).max:
                raise ValueError("deposit region_id is outside the HDF5 int32 range")
            fields = {
                "addition_id": addition_id, "event_id": event_id,
                "particle_id": pid, "process": process,
                "region_id": region_id, "energy_MeV": energy,
                "energy_valid": valid, "kind": kind,
                "start_position_valid": start_valid, "end_position_valid": end_valid,
                "score_position_valid": score_valid,
                "source_record_id": source_record_id,
                "position_semantics": str(deposit.get("position_semantics", "unresolved")),
            }
            for name, value in fields.items():
                dbuf[name].append(value)
            for name, value in (("start_xyz_cm", start), ("end_xyz_cm", end), ("score_xyz_cm", score)):
                dbuf[name].append(value)
            if valid:
                known_total += energy
            else:
                unresolved += 1

        if unresolved and complete:
            raise ValueError("an event containing unresolved deposits cannot be complete")
        if not np.isfinite(known_total):
            raise OverflowError("event deposit total is outside the finite float64 range")
        # Stage every dataset value locally. No buffers, ID sets, or offsets are
        # mutated until all validation above has succeeded.
        staged = {
            "particles": {**pbuf, "track_offsets": track_offsets},
            "tracks": tbuf,
            "deposits": dbuf,
            "events": {
                "event_id": [event_id],
                "particle_offsets": [pstart + len(particles)],
                "deposit_offsets": [dstart + len(deposits)],
                "deposit_total_MeV": [known_total],
                "complete": [bool(complete)],
                "all_absorption": [all_absorption],
                "unresolved_count": [unresolved],
            },
        }
        prior_lengths = {group: {name: len(values) for name, values in fields.items()}
                         for group, fields in self._buffers.items()}
        try:
            for group, fields in staged.items():
                for name, values in fields.items():
                    self._push(group, name, values)
        except Exception:
            for group, fields in self._buffers.items():
                for name, values in fields.items():
                    del values[prior_lengths[group][name]:]
            raise
        self._event_ids.add(event_id)
        self._addition_ids.update(new_addition_ids)
        self._buffered_events += 1
        if self._buffered_events >= 256:
            self._flush()
        return {"particle_start": pstart, "particle_end": self._total("particles", "particle_id"),
                "deposit_start": dstart, "deposit_end": self._total("deposits", "event_id"),
                "deposit_total_MeV": known_total, "unresolved_count": unresolved}

    def close(self, *, commit=True):
        if self.file is None:
            return
        try:
            if self._failed:
                raise RuntimeError("cannot publish an HDF5 writer after a storage failure")
            if commit:
                self._flush()
            self.file.flush()
            self.file.close()
            self.file = None
            if commit:
                # Validate the completed flat arrays before atomically publishing.
                with h5py.File(self.temporary_path, "r") as file:
                    _identify(file, self.temporary_path)
                    if "provenance" not in file or not isinstance(file["provenance"], h5py.Group):
                        raise RuntimeError("missing schema group /provenance")
                os.replace(self.temporary_path, self.path)
            else:
                self.temporary_path.unlink(missing_ok=True)
        except Exception:
            if self.file is not None:
                self.file.close()
                self.file = None
            self.temporary_path.unlink(missing_ok=True)
            raise

    def __enter__(self):
        return self

    def __exit__(self, exc_type, *_):
        self.close(commit=exc_type is None)


def _identify(file, path):
    fmt = file.attrs.get("format", "")
    if _text(fmt) != FORMAT_NAME:
        raise ValueError(f"{path} is not a {FORMAT_NAME} HDF5 file")
    version = file.attrs.get("format_version")
    if version is None:
        raise ValueError(f"{path} has no format_version; legacy hierarchy files use the v1 adapter")
    if int(version) != FORMAT_VERSION:
        raise ValueError(f"unsupported {FORMAT_NAME} format_version={version}")
    if _text(file.attrs.get("schema", "")) != "flat-offset-v2":
        raise ValueError(f"{path} declares v2 with an unknown schema")
    required = {
        "events": ("event_id", "particle_offsets", "deposit_offsets", "deposit_total_MeV",
                   "complete", "all_absorption", "unresolved_count"),
        "particles": ("event_id", "particle_id", "parent_particle_id", "parent_valid", "ityp", "track_offsets"),
        "tracks": ("source_record_id", "ncol", "particle_id", "x_cm", "y_cm", "z_cm", "energy_MeV", "time", "region_id", "time_valid"),
        "deposits": ("addition_id", "event_id", "particle_id", "process", "region_id", "energy_MeV",
                     "energy_valid", "kind", "start_xyz_cm", "end_xyz_cm", "start_position_valid",
                     "end_position_valid", "score_xyz_cm",
                     "score_position_valid", "source_record_id", "position_semantics"),
    }
    for group_name, dataset_names in required.items():
        if group_name not in file or not isinstance(file[group_name], h5py.Group):
            raise ValueError(f"{path} is missing /{group_name}")
        for dataset_name in dataset_names:
            if dataset_name not in file[group_name] or not isinstance(file[group_name][dataset_name], h5py.Dataset):
                raise ValueError(f"{path} is missing /{group_name}/{dataset_name}")

    event_count = len(file["events/event_id"])
    particle_count = len(file["particles/particle_id"])
    track_count = len(file["tracks/ncol"])
    deposit_count = len(file["deposits/event_id"])

    def check_offsets(dataset, expected_length, final_value, label):
        if len(dataset) != expected_length:
            raise ValueError(f"{path} has invalid {label} length")
        previous = 0
        for start in range(0, len(dataset), 1_000_000):
            values = np.asarray(dataset[start:start + 1_000_000], dtype=np.int64)
            if len(values) and (values[0] < previous or np.any(values[1:] < values[:-1])):
                raise ValueError(f"{path} has non-monotonic {label}")
            previous = int(values[-1]) if len(values) else previous
        if int(dataset[0]) != 0 or int(dataset[-1]) != final_value:
            raise ValueError(f"{path} has invalid {label} endpoints")

    check_offsets(file["events/particle_offsets"], event_count + 1, particle_count, "particle_offsets")
    check_offsets(file["events/deposit_offsets"], event_count + 1, deposit_count, "deposit_offsets")
    check_offsets(file["particles/track_offsets"], particle_count + 1, track_count, "track_offsets")
    for group_name, expected in (("events", event_count), ("particles", particle_count),
                                 ("tracks", track_count), ("deposits", deposit_count)):
        group = file[group_name]
        for name, dataset in group.items():
            if name not in ("particle_offsets", "deposit_offsets", "track_offsets") and len(dataset) != expected:
                raise ValueError(f"{path} has inconsistent /{group_name}/{name} length")
    for name in ("start_xyz_cm", "end_xyz_cm", "score_xyz_cm"):
        if file[f"deposits/{name}"].shape != (deposit_count, 3):
            raise ValueError(f"{path} has invalid /deposits/{name} shape")


def iter_deposits(path):
    """Yield ``(event_id, deposit_mapping)`` without collapsing SEGMENT rows."""
    with h5py.File(Path(path), "r") as file:
        _identify(file, path)
        deposits = file["deposits"]
        count = len(deposits["event_id"])
        chunk_size = 65536
        fields = ("addition_id", "event_id", "particle_id", "process", "region_id", "energy_MeV",
                  "energy_valid", "kind", "start_xyz_cm", "end_xyz_cm", "start_position_valid",
                  "end_position_valid", "score_xyz_cm", "score_position_valid", "source_record_id",
                  "position_semantics")
        for base in range(0, count, chunk_size):
            arrays = {name: deposits[name][base:base + chunk_size] for name in fields}
            for row in range(len(arrays["event_id"])):
                kind_id = int(arrays["kind"][row])
                kind = next(name for name, value in DEPOSIT_KINDS.items() if value == kind_id)
                yield int(arrays["event_id"][row]), {
                    "addition_id": _text(arrays["addition_id"][row]),
                    "particle_id": int(arrays["particle_id"][row]),
                    "process": _text(arrays["process"][row]),
                    "region_id": int(arrays["region_id"][row]),
                    "energy_MeV": float(arrays["energy_MeV"][row]) if arrays["energy_valid"][row] else None,
                    "kind": kind,
                    "start_xyz_cm": tuple(float(x) for x in arrays["start_xyz_cm"][row]) if arrays["start_position_valid"][row] else None,
                    "end_xyz_cm": tuple(float(x) for x in arrays["end_xyz_cm"][row]) if arrays["end_position_valid"][row] else None,
                    "score_xyz_cm": tuple(float(x) for x in arrays["score_xyz_cm"][row]) if arrays["score_position_valid"][row] else None,
                    "score_position_valid": bool(arrays["score_position_valid"][row]),
                    "source_record_id": _text(arrays["source_record_id"][row]),
                    "position_semantics": _text(arrays["position_semantics"][row]),
                }


def iter_events(path):
    """Yield event summaries; v2 retains deposit kinds and does not mimic v1."""
    with h5py.File(Path(path), "r") as file:
        _identify(file, path)
        events = file["events"]
        for index, event_id in enumerate(events["event_id"]):
            p0, p1 = map(int, (events["particle_offsets"][index], events["particle_offsets"][index + 1]))
            d0, d1 = map(int, (events["deposit_offsets"][index], events["deposit_offsets"][index + 1]))
            yield int(event_id), {
                "particle_range": (p0, p1), "deposit_range": (d0, d1),
                "deposit_total_MeV": float(events["deposit_total_MeV"][index]),
                "complete": bool(events["complete"][index]),
                "all_absorption": _text(events["all_absorption"][index]),
                "unresolved_count": int(events["unresolved_count"][index]),
            }
