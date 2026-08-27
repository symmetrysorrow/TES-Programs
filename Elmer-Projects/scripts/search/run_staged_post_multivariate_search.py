"""Run a cheap tail screen and send only the best candidates to full search.

The screen and full configurations must use the same variable order, seed, and
sample count so that candidate IDs match.  The screen intentionally uses a
coarser post-pulse grid, a relaxed nonlinear stopping criterion, and a fixed
G0 factor.  Full scoring remains the authoritative result.

Example::

    python scripts/search/run_staged_post_multivariate_search.py --samples 12 --full-top 3
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import post_multivariate_search as search


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_SCREEN_CONFIG = ROOT / "post_multivariate_tail_screen_config.json"
DEFAULT_FULL_CONFIG = ROOT / "post_multivariate_search_config.json"


def absolute(path: Path) -> Path:
    return path if path.is_absolute() else ROOT / path


def score_completed_candidate(
    config: dict, candidate: search.Candidate
) -> bool:
    """Score a pulse that finished before the process failed in scoring."""
    directory = search.candidate_dir(config, candidate.candidate_id)
    metadata_paths = sorted(directory.glob("trial_*/metadata.json"), reverse=True)
    for metadata_path in metadata_paths:
        metadata = search.load_json(metadata_path)
        pulse_cases = metadata.get("pulse_cases", {})
        case_names = {
            str(value["case"])
            for value in pulse_cases.values()
            if isinstance(value, dict) and "case" in value
        }
        if not case_names:
            continue
        if not all(
            all(
                (ROOT / "results" / case_name / search.side_series_file(
                    str(value["series_file"]), side
                )).exists()
                for side in ("L", "R")
            )
            for case_name in case_names
            for value in pulse_cases.values()
            if str(value["case"]) == case_name
        ):
            continue
        score = search.score_candidate_outputs(config, candidate, metadata)
        score.update(
            {
                "project": str((metadata_path.parent / "project.json").relative_to(ROOT)),
                "steady_case": metadata["steady_case"],
                "pulse_cases": pulse_cases,
                "g0_factor": metadata["g0_factor"],
            }
        )
        search.write_json(directory / search.score_filename(config), score)
        search.update_leaderboard(config)
        print(f"[screen] recovered score for {candidate.candidate_id}")
        return True
    return False


def run_screen(
    config: dict,
    candidates: list[search.Candidate],
    *,
    skip_scored: bool,
) -> None:
    score_name = search.score_filename(config)
    selected = candidates
    for index, candidate in enumerate(selected, start=1):
        score_path = search.candidate_dir(config, candidate.candidate_id) / score_name
        if skip_scored and score_path.exists():
            print(f"[screen {index}/{len(selected)}] skip {candidate.candidate_id}")
            continue
        if skip_scored and score_completed_candidate(config, candidate):
            continue
        print(f"[screen {index}/{len(selected)}] {candidate.candidate_id}", flush=True)
        search.run_candidate(config, candidate)


def ranked_screen_candidates(
    config: dict, candidates: list[search.Candidate]
) -> list[search.Candidate]:
    leaderboard = search.output_root(config) / search.leaderboard_filename(config)
    if not leaderboard.exists():
        raise FileNotFoundError(f"screen leaderboard was not created: {leaderboard}")
    by_id = {candidate.candidate_id: candidate for candidate in candidates}
    ranked: list[search.Candidate] = []
    with leaderboard.open("r", newline="", encoding="utf-8") as file:
        import csv

        rows = sorted(csv.DictReader(file), key=lambda row: float(row["objective"]))
    for row in rows:
        candidate = by_id.get(str(row["candidate_id"]))
        if candidate is not None:
            ranked.append(candidate)
    return ranked


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--screen-config", type=Path, default=DEFAULT_SCREEN_CONFIG)
    parser.add_argument("--full-config", type=Path, default=DEFAULT_FULL_CONFIG)
    parser.add_argument("--samples", type=int, default=12)
    parser.add_argument("--seed", type=int, default=20260804)
    parser.add_argument("--full-top", type=int, default=3)
    parser.add_argument("--screen-only", action="store_true")
    parser.add_argument(
        "--full-only",
        action="store_true",
        help="use the existing screen leaderboard without running unfinished screen candidates",
    )
    parser.add_argument("--rerun-screen", action="store_true")
    args = parser.parse_args()
    if args.samples <= 0 or args.full_top <= 0:
        raise SystemExit("--samples and --full-top must be positive")

    screen_config = search.load_json(absolute(args.screen_config))
    full_config = search.load_json(absolute(args.full_config))
    screen_candidates = search.candidates_from_config(
        screen_config, sample_count=args.samples, seed=args.seed
    )
    full_candidates = search.candidates_from_config(
        full_config, sample_count=args.samples, seed=args.seed
    )
    full_by_id = {candidate.candidate_id: candidate for candidate in full_candidates}
    if {candidate.candidate_id for candidate in screen_candidates} != set(full_by_id):
        raise SystemExit(
            "screen and full configurations produce different candidate IDs; "
            "keep variables, seed, and sample count identical"
        )

    if not args.full_only:
        run_screen(
            screen_config,
            screen_candidates,
            skip_scored=not args.rerun_screen,
        )
    ranked = ranked_screen_candidates(screen_config, screen_candidates)
    print("screen ranking:")
    for candidate in ranked[: args.full_top]:
        print(f"  {candidate.candidate_id}: {candidate.factors}")
    if args.screen_only:
        return 0

    selected = ranked[: args.full_top]
    for index, candidate in enumerate(selected, start=1):
        print(f"[full {index}/{len(selected)}] {candidate.candidate_id}", flush=True)
        search.run_candidate(full_config, full_by_id[candidate.candidate_id])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
