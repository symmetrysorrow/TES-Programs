from __future__ import annotations
import argparse, json, sys
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument("--repo-root", required=True, type=Path)
parser.add_argument("--series", required=True, type=Path)
parser.add_argument("--iterations", required=True, type=Path)
parser.add_argument("--audit", required=True, type=Path)
args = parser.parse_args()
sys.path.insert(0, str(args.repo_root / "Elmer-Projects" / "src"))
from support.final_series_flush import flush_final_series_row
result = flush_final_series_row(args.series, args.iterations, True)
args.audit.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
print(json.dumps(result, indent=2))