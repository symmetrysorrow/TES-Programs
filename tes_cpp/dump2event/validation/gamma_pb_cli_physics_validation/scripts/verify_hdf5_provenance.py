#!/usr/bin/env python3
"""Check that CLI JSON and HDF5 v2 provenance hashes agree."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import h5py


HASHES = (
    "physics_contract_sha256",
    "resolved_physics_sha256",
    "physics_config_sha256",
    "phits_input_sha256",
    "source_input_sha256",
    "dumpall_sha256",
    "stopping_table_sha256",
    "phits_executable_sha256",
    "converter_executable_sha256",
    "reconstructor_executable_sha256",
    "transport_exporter_sha256",
    "driver_sha256",
)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--hdf5", type=Path, required=True)
    parser.add_argument("--cli-json", type=Path, required=True)
    args = parser.parse_args()
    raw = args.cli_json.read_bytes()
    encoding = "utf-16" if raw.startswith((b"\xff\xfe", b"\xfe\xff")) else "utf-8-sig"
    result = json.loads(raw.decode(encoding))
    errors: list[str] = []
    with h5py.File(args.hdf5, "r") as h5:
        attrs = h5["provenance"].attrs
        revision = int(attrs.get("provenance_schema_revision", 1))
        for key in HASHES:
            if key not in result:
                if key in attrs:
                    errors.append(f"{key}: HDF5 only")
                continue
            if key not in attrs:
                errors.append(f"{key}: CLI only")
                continue
            actual = attrs[key]
            if isinstance(actual, bytes):
                actual = actual.decode("utf-8")
            if str(actual) != str(result[key]):
                errors.append(f"{key}: CLI/HDF5 mismatch")
        if revision >= 2 and str(attrs.get("physics_config_sha256", "")) != str(
            attrs.get("resolved_physics_sha256", "")
        ):
            errors.append("revision 2 compatibility alias differs from resolved hash")
        resolved_json = attrs.get("resolved_physics_json")
        if isinstance(resolved_json, bytes):
            resolved_json = resolved_json.decode("utf-8")
        event_count = len(h5["events"]["event_id"])
        deposit_count = len(h5["deposits"]["event_id"])
        print(json.dumps({
            "file": str(args.hdf5),
            "revision": revision,
            "events": event_count,
            "deposits": deposit_count,
            "hashes_match": not errors,
            "errors": errors,
            "resolved_physics_json": resolved_json,
        }, indent=2))
    return 1 if errors else 0


if __name__ == "__main__":
    raise SystemExit(main())
