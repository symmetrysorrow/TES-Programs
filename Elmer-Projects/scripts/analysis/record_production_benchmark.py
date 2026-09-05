"""Record the production-size benchmark gate without inventing a run."""
from __future__ import annotations

import argparse
import json
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--candidate-check", type=Path, required=True)
    parser.add_argument("--mesh", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    check = json.loads(args.candidate_check.read_text(encoding="utf-8"))
    header = args.mesh / "mesh.header"
    fields = header.read_text(encoding="utf-8").split() if header.exists() else []
    report = {
        "status": "NOT_RUN",
        "candidate_mesh": str(args.mesh.resolve()),
        "candidate_mesh_counts": {
            "node_count": int(fields[0]) if fields else None,
            "volume_element_count": int(fields[1]) if len(fields) > 1 else None,
        },
        "topology_gate": check.get("status"),
        "reason": "The available production-sized candidate failed the current branch's conformal shared-node interface checker; no CPU/GPU benchmark was run on an invalid route.",
        "required_follow_up": [
            "regenerate the fine production mesh from the current project hash and interface naming",
            "rerun the post-ElmerGrid topology gate",
            "run identical HYPRE CPU/GPU cases and capture setup, linear-solve, and total wall times",
        ],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
