#!/usr/bin/env python3
"""Run the explicitly experimental C++ PHITS 3.341 deposit ledger.

The manifest hashes the exact Dumpall, stopping table, and fixed physics scope.
The executable is deliberately separate from the normal dump2event CLI.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
from datetime import datetime, timezone
from pathlib import Path


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reconstructor", type=Path, required=True)
    parser.add_argument("--dumpall", type=Path, required=True)
    parser.add_argument("--stopping-table", type=Path, required=True)
    parser.add_argument("--table-source", required=True,
                        help="human-readable provenance for the stopping table")
    parser.add_argument("--ledger", type=Path, required=True)
    parser.add_argument("--phits-version", required=True)
    args = parser.parse_args()

    if args.phits_version != "3.341":
        parser.error("this experimental executable is validated only for PHITS 3.341")
    for path in (args.reconstructor, args.dumpall, args.stopping_table):
        if not path.is_file():
            parser.error(f"required file does not exist: {path}")

    config = {
        "phits_version": "3.341",
        "region_id": 1000,
        "region_material": "PB",
        "egs5_medium": 2,
        "electron_cutoff_MeV": 0.001,
        "particle_weights": "unit-only",
        "supported_processes": [
            "NCOL=10/ITYP=12 EGS5 boundary loss",
            "NCOL=13,14/JCOLL=14/ITYP=12 EGS5 continuous loss",
            "NCOL=11/ITYP=12 electron cutoff",
            "dead charged electron secondary above cutoff",
            "strict E == EC zero branch",
        ],
        "status": "experimental-partial; unsupported scoring branches remain unresolved",
    }
    config_bytes = json.dumps(config, sort_keys=True, separators=(",", ":")).encode()
    args.ledger.parent.mkdir(parents=True, exist_ok=True)
    command = [str(args.reconstructor), str(args.dumpall),
               str(args.stopping_table), str(args.ledger)]
    started = datetime.now(timezone.utc)
    completed = subprocess.run(command, check=False, text=True, capture_output=True)
    ended = datetime.now(timezone.utc)
    manifest = {
        "started_utc": started.isoformat(),
        "ended_utc": ended.isoformat(),
        "command": command,
        "exit_code": completed.returncode,
        "stdout": completed.stdout,
        "stderr": completed.stderr,
        "dumpall_sha256": sha256_file(args.dumpall),
        "stopping_table_sha256": sha256_file(args.stopping_table),
        "stopping_table_source": args.table_source,
        "physics_config_sha256": hashlib.sha256(config_bytes).hexdigest(),
        "physics_config": config,
    }
    manifest_path = args.ledger.with_suffix(args.ledger.suffix + ".manifest.json")
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(completed.stdout, end="")
    if completed.stderr:
        print(completed.stderr, end="")
    print(f"manifest={manifest_path}")
    return completed.returncode


if __name__ == "__main__":
    raise SystemExit(main())
