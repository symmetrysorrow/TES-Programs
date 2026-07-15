"""Streaming HDF5 representation of ``dump2event`` event data.

Variable-length particle histories are stored as flat arrays plus offsets.  It
avoids one HDF5 group per event and therefore remains practical for millions
of PHITS histories.
"""

from __future__ import annotations

import json
from pathlib import Path

import h5py
import numpy as np

FORMAT_NAME = "tes-dump2event"
FORMAT_VERSION = 1


def iter_json_object_items(path, chunk_size=1024 * 1024):
    """Stream top-level JSON object members without retaining the document."""
    decoder = json.JSONDecoder()
    with open(path, encoding="utf-8") as file:
        buffer = ""

        def read_more():
            nonlocal buffer
            chunk = file.read(chunk_size)
            if not chunk:
                raise ValueError(f"unexpected end of JSON input: {path}")
            buffer += chunk

        def parse_value():
            nonlocal buffer
            while True:
                buffer = buffer.lstrip()
                try:
                    value, end = decoder.raw_decode(buffer)
                except json.JSONDecodeError:
                    read_more()
                    continue
                buffer = buffer[end:]
                return value

        while not buffer:
            read_more()
        buffer = buffer.lstrip()
        if not buffer.startswith("{"):
            raise ValueError(f"expected a JSON object at the root of: {path}")
        buffer = buffer[1:]
        while True:
            while not buffer.lstrip():
                read_more()
            buffer = buffer.lstrip()
            if buffer.startswith("}"):
                return
            if buffer.startswith(","):
                buffer = buffer[1:]
                continue
            key = parse_value()
            while not buffer.lstrip():
                read_more()
            buffer = buffer.lstrip()
            if not buffer.startswith(":"):
                raise ValueError(f"expected ':' after key {key!r}")
            buffer = buffer[1:]
            yield key, parse_value()


class EventWriter:
    def __init__(self, path):
        self.file = h5py.File(Path(path), "w")
        self.file.attrs["format"] = FORMAT_NAME
        self.file.attrs["format_version"] = FORMAT_VERSION
        self.events = self._dataset("event_id", np.int64)
        self.event_offsets = self._dataset("event_particle_offsets", np.int64)
        self.event_offsets.resize(1, axis=0)
        self.event_offsets[0] = 0
        self.particle_ids = self._dataset("particle_id", np.int64)
        self.particle_types = self._dataset("ityp", np.int32)
        self.collision_offsets = self._dataset("particle_collision_offsets", np.int64)
        self.deposit_offsets = self._dataset("particle_deposit_offsets", np.int64)
        for offsets in (self.collision_offsets, self.deposit_offsets):
            offsets.resize(1, axis=0)
            offsets[0] = 0
        self.collision = {name: self._dataset(name, np.float64) for name in ("x", "y", "z", "E")}
        self.deposit = {name: self._dataset(name, np.float64) for name in ("x_deposit", "y_deposit", "z_deposit", "E_deposit")}

    def _dataset(self, name, dtype):
        return self.file.create_dataset(name, shape=(0,), maxshape=(None,), dtype=dtype,
                                        chunks=True, compression="gzip")

    @staticmethod
    def _append(dataset, values):
        values = np.asarray(values, dtype=dataset.dtype)
        start = len(dataset)
        dataset.resize(start + len(values), axis=0)
        dataset[start:] = values

    def append(self, event_id, particles):
        event_index = len(self.events)
        self._append(self.events, [int(event_id)])
        for particle_id, particle in particles.items():
            collision_size = len(particle.get("x", []))
            deposit_size = len(particle.get("x_deposit", []))
            for field in ("y", "z", "E"):
                if len(particle.get(field, [])) != collision_size:
                    raise ValueError(f"event {event_id}, particle {particle_id}: collision arrays have different lengths")
            for field in ("y_deposit", "z_deposit", "E_deposit"):
                if len(particle.get(field, [])) != deposit_size:
                    raise ValueError(f"event {event_id}, particle {particle_id}: deposit arrays have different lengths")
            self._append(self.particle_ids, [int(particle_id)])
            self._append(self.particle_types, [int(particle.get("ityp", 0))])
            self._append(self.collision_offsets, [len(self.collision["x"]) + collision_size])
            self._append(self.deposit_offsets, [len(self.deposit["x_deposit"]) + deposit_size])
            for field, dataset in self.collision.items():
                self._append(dataset, particle.get(field, []))
            for field, dataset in self.deposit.items():
                self._append(dataset, particle.get(field, []))
        self._append(self.event_offsets, [len(self.particle_ids)])

    def close(self):
        self.file.close()

    def __enter__(self):
        return self

    def __exit__(self, *_):
        self.close()


def json_to_hdf5(json_path, hdf5_path):
    with EventWriter(hdf5_path) as writer:
        for event_id, particles in iter_json_object_items(json_path):
            if not isinstance(particles, dict):
                raise ValueError(f"event {event_id} is not an object")
            writer.append(event_id, particles)


def iter_events(path):
    """Yield ``(event_id, particles)`` in the legacy event.json shape."""
    with h5py.File(path, "r") as file:
        if "events" in file:
            for event_id, event in file["events"].items():
                particles = {}
                for particle_id, particle in event.items():
                    particles[particle_id] = {"ityp": int(particle.attrs["ityp"]),
                        **{name: particle[name][:].tolist() for name in ("x", "y", "z", "E", "x_deposit", "y_deposit", "z_deposit", "E_deposit")}}
                yield event_id, particles
            return
        if file.attrs.get("format", "") != FORMAT_NAME:
            raise ValueError(f"{path} is not a {FORMAT_NAME} HDF5 file")
        ids = file["event_id"]
        event_offsets = file["event_particle_offsets"]
        particle_ids = file["particle_id"]
        particle_types = file["ityp"]
        collision_offsets = file["particle_collision_offsets"]
        deposit_offsets = file["particle_deposit_offsets"]
        for event_index, event_id in enumerate(ids):
            particles = {}
            for particle_index in range(event_offsets[event_index], event_offsets[event_index + 1]):
                c0, c1 = collision_offsets[particle_index], collision_offsets[particle_index + 1]
                d0, d1 = deposit_offsets[particle_index], deposit_offsets[particle_index + 1]
                particles[str(particle_ids[particle_index])] = {
                    "ityp": int(particle_types[particle_index]),
                    **{field: file[field][c0:c1].tolist() for field in ("x", "y", "z", "E")},
                    **{field: file[field][d0:d1].tolist() for field in ("x_deposit", "y_deposit", "z_deposit", "E_deposit")},
                }
            yield str(event_id), particles
