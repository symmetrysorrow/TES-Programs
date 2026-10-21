"""Strict, additive-only comparison of reconstructed PHITS tally increments.

This tool validates an explicitly matched tally-increment ledger. It does NOT
establish that a tally scoring coordinate is a microscopic energy-loss point,
or reconstruct any physics from Dumpall by itself.

No PHITS data, stopping-power table, or proprietary Fortran source is bundled.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
import sys
from typing import Dict, List, Optional, Tuple


class InvalidLedger(ValueError):
    """Input is malformed, ambiguous, or unsafe to compare."""


def _finite(value: str, where: str) -> float:
    try:
        result = float(value)
    except (ValueError, TypeError) as exc:
        raise InvalidLedger(f"{where}: invalid number {value!r}") from exc
    if not math.isfinite(result):
        raise InvalidLedger(f"{where}: non-finite number {value!r}")
    return result


def load_ledger(path: Path, *, key_column: str, energy_column: str,
                history_column: Optional[str], region_column: Optional[str],
                position_columns: Optional[Tuple[str, str, str]]) -> Dict[str, dict]:
    required = {key_column, energy_column}
    if history_column:
        required.add(history_column)
    if region_column:
        required.add(region_column)
    if position_columns:
        required.update(position_columns)
    result: Dict[str, dict] = {}
    try:
        with path.open('r', encoding='utf-8-sig', newline='') as stream:
            reader = csv.DictReader(stream)
            if not reader.fieldnames:
                raise InvalidLedger(f"{path}: missing CSV header")
            missing = required.difference(reader.fieldnames)
            if missing:
                raise InvalidLedger(f"{path}: missing columns: {', '.join(sorted(missing))}")
            for line, row in enumerate(reader, start=2):
                if None in row or any(value is None for value in row.values()):
                    raise InvalidLedger(f"{path}:{line}: ragged CSV row")
                key = row[key_column].strip()
                if not key:
                    raise InvalidLedger(f"{path}:{line}: blank tally-increment key")
                if key in result:
                    raise InvalidLedger(f"{path}:{line}: duplicate tally-increment key {key!r}")
                energy = _finite(row[energy_column], f"{path}:{line} {energy_column}")
                if energy < 0:
                    raise InvalidLedger(f"{path}:{line}: negative physical tally contribution {energy}")
                history = row[history_column].strip() if history_column else None
                region = row[region_column].strip() if region_column else None
                if history_column and not history:
                    raise InvalidLedger(f"{path}:{line}: blank history identifier")
                if region_column and not region:
                    raise InvalidLedger(f"{path}:{line}: blank region identifier")
                xyz = tuple(_finite(row[c], f"{path}:{line} {c}") for c in position_columns) if position_columns else None
                result[key] = {'energy': energy, 'history': history, 'region': region, 'xyz': xyz}
    except OSError as exc:
        raise InvalidLedger(f"could not read {path}: {exc}") from exc
    if not result:
        raise InvalidLedger(f"{path}: empty ledger; no comparison performed")
    return result


def compare(reference: Dict[str, dict], candidate: Dict[str, dict], *,
            energy_abs_tol: float, position_abs_tol: float) -> dict:
    keys_ref, keys_candidate = set(reference), set(candidate)
    missing, unexpected = sorted(keys_ref - keys_candidate), sorted(keys_candidate - keys_ref)
    common = sorted(keys_ref & keys_candidate)
    failures: List[dict] = []
    max_energy_error = 0.0
    max_position_error = 0.0
    ref_histories: Dict[str, List[float]] = {}
    cand_histories: Dict[str, List[float]] = {}
    for key in common:
        ref, got = reference[key], candidate[key]
        delta = abs(ref['energy'] - got['energy'])
        max_energy_error = max(max_energy_error, delta)
        problems = []
        if delta > energy_abs_tol:
            problems.append('energy')
        if ref['history'] != got['history']:
            problems.append('history')
        if ref['region'] != got['region']:
            problems.append('region')
        if ref['xyz'] is not None:
            coord_err = max(abs(a - b) for a, b in zip(ref['xyz'], got['xyz']))
            max_position_error = max(max_position_error, coord_err)
            if coord_err > position_abs_tol:
                problems.append('position')
        if problems and len(failures) < 20:
            failures.append({'key': key, 'checks': problems,
                             'reference_energy': ref['energy'], 'candidate_energy': got['energy']})
        if ref['history'] is not None:
            ref_histories.setdefault(ref['history'], []).append(ref['energy'])
            cand_histories.setdefault(got['history'], []).append(got['energy'])
    # Include unmatched rows in the per-history ledger: otherwise missing
    # increments can be hidden by comparing only the intersection.
    if reference and next(iter(reference.values()))['history'] is not None:
        ref_histories.clear()
        cand_histories.clear()
        for item in reference.values():
            ref_histories.setdefault(item['history'], []).append(item['energy'])
        for item in candidate.values():
            cand_histories.setdefault(item['history'], []).append(item['energy'])
    ref_per_history = {k: math.fsum(v) for k, v in ref_histories.items()}
    cand_per_history = {k: math.fsum(v) for k, v in cand_histories.items()}
    all_histories = set(ref_per_history) | set(cand_per_history)
    hist_diff = max((abs(ref_per_history.get(h, 0.0) - cand_per_history.get(h, 0.0))
                     for h in all_histories), default=0.0)
    hist_mismatch = sum(abs(ref_per_history.get(h, 0.0) - cand_per_history.get(h, 0.0)) > energy_abs_tol
                        for h in all_histories)
    # Do not make a successful total sum mask a mismatched constituent entry.
    mismatched_entries = sum(
        (abs(reference[k]['energy'] - candidate[k]['energy']) > energy_abs_tol
         or reference[k]['history'] != candidate[k]['history']
         or reference[k]['region'] != candidate[k]['region']
         or (reference[k]['xyz'] is not None and
             max(abs(a - b) for a, b in zip(reference[k]['xyz'], candidate[k]['xyz'])) > position_abs_tol))
        for k in common
    )
    passed = not missing and not unexpected and not mismatched_entries and not hist_mismatch
    return {
        'passed': passed, 'reference_count': len(reference), 'candidate_count': len(candidate),
        'matched_count': len(common), 'missing_count': len(missing),
        'unexpected_count': len(unexpected), 'mismatched_count': mismatched_entries,
        'reference_zero_count': sum(v['energy'] == 0 for v in reference.values()),
        'candidate_zero_count': sum(v['energy'] == 0 for v in candidate.values()),
        'reference_total_mev': math.fsum(v['energy'] for v in reference.values()),
        'candidate_total_mev': math.fsum(v['energy'] for v in candidate.values()),
        'max_energy_abs_error_mev': max_energy_error,
        'max_position_abs_error': max_position_error,
        'history_count': len(all_histories), 'mismatched_history_count': hist_mismatch,
        'max_history_abs_error_mev': hist_diff,
        'missing_example_keys': missing[:20], 'unexpected_example_keys': unexpected[:20],
        'mismatch_examples': failures,
    }


def main(argv: Optional[List[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('reference_csv', type=Path)
    parser.add_argument('candidate_csv', type=Path)
    parser.add_argument('--key-column', required=True, help='stable, unique ID for each PHITS tally addition')
    parser.add_argument('--reference-energy-column', required=True)
    parser.add_argument('--candidate-energy-column', required=True)
    parser.add_argument('--history-column', help='same history-ID column in both CSVs')
    parser.add_argument('--region-column', help='same scoring-region column in both CSVs')
    parser.add_argument('--position-columns', nargs=3, metavar=('X', 'Y', 'Z'),
                        help='scoring coordinates, not claimed physical microscopic positions')
    parser.add_argument('--energy-abs-tol', type=float, default=1e-12)
    parser.add_argument('--position-abs-tol', type=float, default=1e-10)
    parser.add_argument('--json-report', type=Path)
    args = parser.parse_args(argv)
    if not (math.isfinite(args.energy_abs_tol) and args.energy_abs_tol >= 0 and
            math.isfinite(args.position_abs_tol) and args.position_abs_tol >= 0):
        parser.error('tolerances must be finite and nonnegative')
    try:
        opts = dict(key_column=args.key_column, history_column=args.history_column,
                    region_column=args.region_column, position_columns=args.position_columns)
        ref = load_ledger(args.reference_csv, energy_column=args.reference_energy_column, **opts)
        cand = load_ledger(args.candidate_csv, energy_column=args.candidate_energy_column, **opts)
        result = compare(ref, cand, energy_abs_tol=args.energy_abs_tol,
                         position_abs_tol=args.position_abs_tol)
    except InvalidLedger as exc:
        print(f'INVALID: {exc}', file=sys.stderr)
        return 2
    output = json.dumps(result, indent=2, ensure_ascii=False, allow_nan=False) + '\n'
    if args.json_report:
        try:
            args.json_report.write_text(output, encoding='utf-8')
        except OSError as exc:
            print(f'INVALID: cannot write report: {exc}', file=sys.stderr)
            return 2
    print(output, end='')
    return 0 if result['passed'] else 1


if __name__ == '__main__':
    sys.exit(main())
