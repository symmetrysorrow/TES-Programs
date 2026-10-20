"""Appendable flat HDF5 v2 schema for PHITS transport and deposit records.

This module deliberately keeps transport segments distinct from point scores.
It is a storage/API layer; it does not reconstruct physics or assign TES blocks.
"""

from __future__ import annotations

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
        for key in REQUIRED_PROVENANCE:
            self.provenance.attrs[key] = str(provenance[key])
        for key in ("dumpall_sha256", "physics_config_sha256", "stopping_table_sha256"):
            if not re.fullmatch(r"[0-9a-fA-F]{64}", str(provenance[key])):
                self.file.close()
                self.temporary_path.unlink(missing_ok=True)
                raise ValueError(f"{key} must be a 64-character SHA-256 hex digest")

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

    def _total(self, group: str, name: str) -> int:
        return len(getattr(self, group)[name]) + len(self._buffers[group][name])

    def _push(self, group: str, name: str, values) -> None:
        self._buffers[group][name].extend(values)

    def _flush(self) -> None:
        for group_name, fields in self._buffers.items():
            group = getattr(self, group_name)
            for name, values in fields.items():
                if values:
                    _append(group[name], values)
                    fields[name] = []
        self._buffered_events = 0
        self.file.flush()

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
        if all_absorption not in ABSORPTION_STATES:
            raise ValueError(f"invalid all_absorption state: {all_absorption}")
        if int(event_id) in self._event_ids:
            raise ValueError(f"duplicate event_id {event_id}")
        pstart = self._total("particles", "particle_id")
        dstart = self._total("deposits", "event_id")
        particle_ids = set()
        known_total = 0.0
        unresolved = 0
        pbuf = {name: [] for name in ("event_id", "particle_id", "parent_particle_id", "parent_valid", "ityp")}
        tbuf = {name: [] for name in ("source_record_id", "ncol", "particle_id", "x_cm", "y_cm", "z_cm",
                                      "energy_MeV", "time", "region_id", "time_valid")}
        dbuf = {name: [] for name in ("addition_id", "event_id", "particle_id", "process", "region_id",
                                      "energy_MeV", "energy_valid", "kind", "start_xyz_cm", "end_xyz_cm",
                                      "score_xyz_cm", "score_position_valid", "source_record_id", "position_semantics")}
        track_offsets = []
        for particle in particles:
            pid = int(particle["particle_id"])
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
                if not np.isfinite(energy):
                    raise ValueError("track energy must be finite")
                tbuf["source_record_id"].append(str(track["source_record_id"]))
                tbuf["ncol"].append(int(track["ncol"]))
                tbuf["particle_id"].append(pid)
                for axis, value in zip(("x_cm", "y_cm", "z_cm"), xyz):
                    tbuf[axis].append(value)
                time_valid = track.get("time") is not None
                time = float(track["time"]) if time_valid else 0.0
                if not np.isfinite(time):
                    raise ValueError("track time must be finite")
                tbuf["time"].append(time)
                tbuf["time_valid"].append(time_valid)
                tbuf["region_id"].append(int(track["region_id"]))
            track_offsets.append(self._total("tracks", "ncol") + len(tbuf["ncol"]))

        for name, values in pbuf.items():
            self._push("particles", name, values)
        for name, values in tbuf.items():
            self._push("tracks", name, values)
        self._push("particles", "track_offsets", track_offsets)

        for deposit in deposits:
            addition_id = str(deposit["addition_id"])
            if addition_id in self._addition_ids:
                raise ValueError(f"duplicate addition_id {addition_id}")
            self._addition_ids.add(addition_id)
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
            score_valid = deposit.get("score_xyz_cm") is not None
            score = _vector(deposit.get("score_xyz_cm"), "score_xyz_cm")
            pid = int(deposit.get("particle_id", -1))
            if pid >= 0 and pid not in particle_ids:
                raise ValueError(f"deposit references absent particle_id {pid} in event {event_id}")
            fields = {
                "addition_id": str(deposit["addition_id"]), "event_id": event_id,
                "particle_id": pid, "process": str(deposit["process"]),
                "region_id": int(deposit["region_id"]), "energy_MeV": energy,
                "energy_valid": valid, "kind": kind,
                "score_position_valid": score_valid,
                "source_record_id": str(deposit["source_record_id"]),
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

        for name, values in dbuf.items():
            self._push("deposits", name, values)

        if unresolved and complete:
            raise ValueError("an event containing unresolved deposits cannot be complete")
        self._event_ids.add(int(event_id))
        self._push("events", "event_id", [event_id])
        self._push("events", "particle_offsets", [self._total("particles", "particle_id")])
        self._push("events", "deposit_offsets", [self._total("deposits", "event_id")])
        self._push("events", "deposit_total_MeV", [known_total])
        self._push("events", "complete", [bool(complete)])
        self._push("events", "all_absorption", [all_absorption])
        self._push("events", "unresolved_count", [unresolved])
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
            if commit:
                self._flush()
            self.file.flush()
            self.file.close()
            self.file = None
            if commit:
                # Validate the completed flat arrays before atomically publishing.
                with h5py.File(self.temporary_path, "r") as file:
                    n_events = len(file["events/event_id"])
                    if len(file["events/particle_offsets"]) != n_events + 1 or len(file["events/deposit_offsets"]) != n_events + 1:
                        raise RuntimeError("invalid event offset arrays")
                    if len(file["particles/track_offsets"]) != len(file["particles/particle_id"]) + 1:
                        raise RuntimeError("invalid track offset array")
                    for group in ("events", "particles", "tracks", "deposits", "provenance"):
                        if group not in file:
                            raise RuntimeError(f"missing schema group /{group}")
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


def iter_deposits(path):
    """Yield ``(event_id, deposit_mapping)`` without collapsing SEGMENT rows."""
    with h5py.File(Path(path), "r") as file:
        _identify(file, path)
        events = file["events"]
        deposits = file["deposits"]
        for index, event_id in enumerate(events["event_id"]):
            start, end = map(int, (events["deposit_offsets"][index], events["deposit_offsets"][index + 1]))
            for row in range(start, end):
                kind_id = int(deposits["kind"][row])
                kind = next(name for name, value in DEPOSIT_KINDS.items() if value == kind_id)
                yield int(event_id), {
                    "addition_id": _text(deposits["addition_id"][row]),
                    "particle_id": int(deposits["particle_id"][row]),
                    "process": _text(deposits["process"][row]),
                    "region_id": int(deposits["region_id"][row]),
                    "energy_MeV": float(deposits["energy_MeV"][row]) if deposits["energy_valid"][row] else None,
                    "kind": kind,
                    "start_xyz_cm": tuple(float(x) for x in deposits["start_xyz_cm"][row]),
                    "end_xyz_cm": tuple(float(x) for x in deposits["end_xyz_cm"][row]),
                    "score_xyz_cm": tuple(float(x) for x in deposits["score_xyz_cm"][row]) if deposits["score_position_valid"][row] else None,
                    "score_position_valid": bool(deposits["score_position_valid"][row]),
                    "source_record_id": _text(deposits["source_record_id"][row]),
                    "position_semantics": _text(deposits["position_semantics"][row]),
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
