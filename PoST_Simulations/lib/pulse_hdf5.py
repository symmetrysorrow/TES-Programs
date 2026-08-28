"""Compact, streaming HDF5 storage for TES pulse waveforms.

Schema (format version 1): ``time`` is a shared one-dimensional dataset and
``event_id``, ``ch0`` and ``ch1`` are aligned along their first dimension.
Waveforms are appended one at a time, so generating large event sets does not
require retaining them in memory.
"""

from __future__ import annotations

import json
from pathlib import Path

import h5py
import numpy as np

FORMAT_NAME = "tes-pulses"
FORMAT_VERSION = 1


class PulseWriter:
    def __init__(self, path, time, input_parameters=None):
        self.path = Path(path)
        time = np.asarray(time, dtype=np.float64)
        if time.ndim != 1 or len(time) == 0:
            raise ValueError("time must be a non-empty one-dimensional array")
        self.samples = len(time)
        self.file = h5py.File(self.path, "w")
        self.file.attrs["format"] = FORMAT_NAME
        self.file.attrs["format_version"] = FORMAT_VERSION
        if input_parameters is not None:
            self.file.attrs["input_json"] = json.dumps(input_parameters, separators=(",", ":"))
        self.file.create_dataset("time", data=time)
        self.ids = self.file.create_dataset(
            "event_id", shape=(0,), maxshape=(None,), dtype=h5py.string_dtype("utf-8")
        )
        options = dict(shape=(0, self.samples), maxshape=(None, self.samples),
                       dtype=np.float64, chunks=(1, self.samples), compression="gzip")
        self.ch0 = self.file.create_dataset("ch0", **options)
        self.ch1 = self.file.create_dataset("ch1", **options)

    def append(self, event_id, ch0, ch1):
        ch0 = np.asarray(ch0, dtype=np.float64)
        ch1 = np.asarray(ch1, dtype=np.float64)
        if ch0.shape != (self.samples,) or ch1.shape != (self.samples,):
            raise ValueError(f"pulse {event_id!r} does not match {self.samples} samples")
        index = len(self.ids)
        for dataset in (self.ids, self.ch0, self.ch1):
            dataset.resize(index + 1, axis=0)
        self.ids[index] = str(event_id)
        self.ch0[index] = ch0
        self.ch1[index] = ch1

    def close(self):
        self.file.close()

    def __enter__(self):
        return self

    def __exit__(self, *_):
        self.close()


def iter_pulse_items(path, event_ids=None):
    """Yield one ``(event_id, {ch0, ch1})`` pair at a time.

    When *event_ids* is supplied, skip non-selected rows before reading the
    two large waveform datasets.  This is important for FullEnergy analyses,
    where only a small subset of the stored events is needed.
    """
    selected_ids = None if event_ids is None else {str(event_id) for event_id in event_ids}
    with h5py.File(path, "r") as file:
        if file.attrs.get("format", "") != FORMAT_NAME:
            raise ValueError(f"{path} is not a {FORMAT_NAME} HDF5 file")
        ids, ch0, ch1 = file["event_id"], file["ch0"], file["ch1"]
        if not (len(ids) == len(ch0) == len(ch1)):
            raise ValueError(f"{path} has inconsistent pulse dataset lengths")
        for index, event_id in enumerate(ids.asstr()):
            if selected_ids is not None and event_id not in selected_ids:
                continue
            yield event_id, {"ch0": ch0[index], "ch1": ch1[index]}


def read_pulse(path, event_id):
    with h5py.File(path, "r") as file:
        ids = file["event_id"].asstr()[:]
        matches = np.flatnonzero(ids == str(event_id))
        if not len(matches):
            raise KeyError(f"Event ID {event_id!r} was not found in {path}")
        index = int(matches[0])
        return {"ch0": file["ch0"][index], "ch1": file["ch1"][index]}


def read_time(path):
    with h5py.File(path, "r") as file:
        return file["time"][:]


def read_all_pulses(path):
    time = read_time(path)
    return {
        int(event_id): {"position": int(event_id), "time": time, **pulse}
        for event_id, pulse in iter_pulse_items(path)
    }
