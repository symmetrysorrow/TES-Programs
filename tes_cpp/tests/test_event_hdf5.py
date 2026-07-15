"""Smoke test for streaming dump2event HDF5 conversion."""

from __future__ import annotations

import json
import tempfile
from pathlib import Path

from tes_cpp.event_hdf5 import iter_events, json_to_hdf5


def main() -> None:
    document = {
        "11": {
            "1": {"ityp": 14, "x": [1.0], "y": [2.0], "z": [3.0], "E": [4.0],
                  "x_deposit": [1.5], "y_deposit": [2.5], "z_deposit": [3.5], "E_deposit": [0.4]},
            "2": {"ityp": 12, "x": [], "y": [], "z": [], "E": [],
                  "x_deposit": [], "y_deposit": [], "z_deposit": [], "E_deposit": []},
        },
        "12": {},
    }
    with tempfile.TemporaryDirectory() as directory:
        json_path = Path(directory) / "event.json"
        hdf5_path = Path(directory) / "event.h5"
        json_path.write_text(json.dumps(document), encoding="utf-8")
        json_to_hdf5(json_path, hdf5_path)
        assert dict(iter_events(hdf5_path)) == document
    print("PASS: dump2event HDF5 round-trip")


if __name__ == "__main__":
    main()
