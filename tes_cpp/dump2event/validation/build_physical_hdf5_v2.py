#!/usr/bin/env python3
"""Stream the validated Pb-only experimental ledger and raw tracks into HDF5 v2.

This is an opt-in research utility. It does not call the ordinary converter and
does not claim support for physical processes outside the existing Pb ledger.
"""
from __future__ import annotations

import argparse
import csv
from decimal import Decimal, InvalidOperation
import hashlib
import itertools
import json
import math
import os
import re
import subprocess
import sys
import time
import tempfile
from pathlib import Path

# Permit the installed/standalone validation driver to import the repository's
# small HDF5 v2 writer without requiring callers to edit PYTHONPATH.
_python_root = Path(__file__).resolve().parents[2] / "python"
if _python_root.is_dir() and str(_python_root) not in sys.path:
    sys.path.insert(0, str(_python_root))

from tes_cpp.event_hdf5_v2 import EventV2Writer, iter_deposits, iter_events


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def canonical_json(value: object) -> bytes:
    """Stable UTF-8 representation for semantic configuration hashes."""
    return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=False).encode("utf-8")


MAX_SOURCE_GAMMA_ENERGY_MEV = 2.0
PHITS_3341_EGS5_EMIN_DEFAULTS_MEV = {12: 0.1, 13: 0.1, 14: 0.001}
VALIDATED_EGS5_CUTOFFS_MEV = {12: 0.001, 13: 0.001, 14: 0.001}


def source_definition(source_lines: list[str]) -> dict:
    """Accept only unit-weight discrete photon lines in (0, 2] MeV.

    This is a source syntax gate, NOT an independent PHITS scoring validation
    at every energy. The verified Pb HATCH material/cutoff and per-process
    reconstruction guards still apply separately.
    """
    values: dict[str, str] = {}
    energies: list[float] = []
    for line in source_lines:
        if "=" in line:
            key, value = line.split("=", 1)
            key = re.sub(r"\s+", "", key).lower()
            if key in values:
                raise ValueError(f"duplicate photon source setting: {key}")
            values[key] = value.strip()
            continue
        fields = line.split()
        if len(fields) != 2:
            raise ValueError(f"unsupported photon source energy row: {line!r}")
        try:
            energy, intensity = map(float, fields)
        except ValueError as error:
            raise ValueError(f"invalid photon source energy row: {line!r}") from error
        if not (math.isfinite(energy) and 0.0 < energy <= MAX_SOURCE_GAMMA_ENERGY_MEV):
            raise ValueError("photon source energy must be finite, >0 and <=2.0 MeV")
        if not (math.isfinite(intensity) and intensity == 1.0):
            raise ValueError("only unit-intensity photon lines are supported")
        energies.append(energy)
    if not energies:
        raise ValueError("photon source must declare at least one discrete energy line")
    # Omitted totfact is not inferred as a PHITS tally normalization; it is not
    # consulted by the Dumpall score reconstruction. If explicitly supplied,
    # retain the unit-intensity boundary of the validated implementation.
    try:
        totfact = float(values.get("totfact", "1"))
        wgt = float(values.get("wgt", "1"))
    except ValueError as error:
        raise ValueError("photon source totfact and wgt must be numeric unit weights") from error
    if not (math.isfinite(totfact) and totfact == 1.0 and
            math.isfinite(wgt) and wgt == 1.0):
        raise ValueError("photon source requires unit totfact=1 and wgt=1")

    # Preserve historical Co-60 inputs and provenance names. Earlier versions
    # did not require these three fields for this exact two-line source, so
    # omitted fields are accepted ONLY for the original validated Co-60 pair.
    co60 = energies == [1.173, 1.332]
    # ne is metadata for the line count, not an input to reconstruction.
    # If PHITS accepted a source without explicit ne, the line list already
    # defines the count; if supplied, reject contradictory values.
    required = ("proj", "e-type")
    if not co60 and any(key not in values for key in required):
        raise ValueError("photon source requires proj=photon and e-type=8")
    if "proj" in values and values["proj"].lower() != "photon":
        raise ValueError("photon source requires proj=photon")
    if "e-type" in values and values["e-type"] != "8":
        raise ValueError("photon source requires e-type=8")
    if "ne" in values:
        try:
            declared_ne = int(values["ne"])
        except ValueError as error:
            raise ValueError("photon source ne must be an integer") from error
        if declared_ne != len(energies):
            raise ValueError("photon source ne must equal number of energy lines")

    if co60:
        mode = "co60-lines"
    elif energies == [0.511]:
        mode = "mono-gamma-511keV"
    elif energies == [2.0]:
        mode = "mono-gamma-2MeV"
    else:
        mode = "mono-gamma" if len(energies) == 1 else "discrete-gamma-lines"
    return {"source_mode": mode, "source_particle": "photon",
            "source_energies_MeV": energies}

def resolved_esmin_MeV(source_input: Path, phits_version: str) -> float:
    """Resolve PHITS's global [Parameters] esmin, using only the pinned 3.341 default."""
    if phits_version != "3.341":
        raise ValueError("dead-secondary score gate is validated only for PHITS 3.341")
    active = ""
    params: dict[str, str] = {}
    for raw in source_input.read_text(encoding="utf-8-sig", errors="replace").splitlines():
        line = raw.split("#", 1)[0].strip()
        if not line:
            continue
        if line.startswith("[") and line.endswith("]"):
            active = re.sub(r"\s+", "", line[1:-1]).lower()
        elif active == "parameters" and "=" in line:
            key, value = line.split("=", 1)
            params[re.sub(r"\s+", "", key).lower()] = value.strip()
    explicit = params.get("esmin")
    if explicit is None:
        # PHITS 3.341 read02.f maps the `esmin` input name to parz(162).
        # The checked-in default table initializes parz(162) to exactly 1.0e-3 MeV (confirmed from read02.f data parz(161:170), element 2).
        return 1.0e-3
    try:
        value = float(explicit.split()[0])
    except (ValueError, IndexError) as error:
        raise ValueError(f"invalid [Parameters] esmin: {explicit!r}") from error
    if not (value > 0.0 and value < float("inf")):
        raise ValueError("[Parameters] esmin must be finite and positive")
    return value


def _parameter_number(raw: str, name: str) -> float:
    """Parse a PHITS numeric token, including Fortran D exponents."""
    try:
        token = raw.split()[0].replace("D", "E").replace("d", "e")
        value = float(token)
    except (ValueError, IndexError) as error:
        raise ValueError(f"invalid [Parameters] {name}: {raw!r}") from error
    if not math.isfinite(value):
        raise ValueError(f"[Parameters] {name} must be finite")
    return value


def resolved_egs5_cutoffs(source_input: Path, phits_version: str) -> dict:
    """Resolve EGS5 emin(12:14) from explicit values or PHITS 3.341 defaults."""
    if phits_version != "3.341":
        raise ValueError("EGS5 cutoff defaults are pinned only for PHITS 3.341")
    params = _input_values(source_input, "parameters")
    raw_negs = params.get("negs")
    negs = _parameter_number(raw_negs, "negs") if raw_negs is not None else -1.0
    if negs not in (1.0, 2.0):
        source = "explicit" if raw_negs is not None else "PHITS default"
        raise ValueError(f"EGS5 cutoff resolution requires negs=1 or 2; effective value is {negs:g} ({source})")
    values: dict[int, float] = {}
    sources: dict[str, str] = {}
    for index, default in PHITS_3341_EGS5_EMIN_DEFAULTS_MEV.items():
        key = f"emin({index})"
        raw = params.get(key)
        values[index] = default if raw is None else _parameter_number(raw, key)
        if values[index] <= 0.0:
            raise ValueError(f"[Parameters] {key} must be positive; found {values[index]!r}")
        sources[key] = "phits_default" if raw is None else "explicit"
    # PHITS 3.341 enforces emin(12) == emin(13). If they differ, the run
    # reports "emin(12) is adjusted to emin(13)" and uses the positron value
    # for both charged cutoffs. Resolve the effective value, not the raw
    # independent default for emin(12).
    if values[12] != values[13]:
        values[12] = values[13]
        sources["emin(12)"] = "phits_adjusted_to_emin13"
    return {
        "electron_cutoff_MeV": values[12],
        "positron_cutoff_MeV": values[13],
        "photon_cutoff_MeV": values[14],
        "emin_sources": sources,
    }


def _resolved_run_path(input_path: Path, raw_path: str) -> Path:
    reference = Path(raw_path.strip().strip("\"'"))
    return reference.resolve() if reference.is_absolute() else (input_path.parent / reference).resolve()


def _log_parameter(log_text: str, key: str) -> str | None:
    pattern = rf"(?im)^\s*{re.escape(key)}\s*=\s*(.*?)\s*(?:#.*)?$"
    matches = re.findall(pattern, log_text)
    return matches[-1].strip().strip("\"'") if matches else None


def material_settings(contract: dict) -> dict:
    """Resolve explicit Pb/Bi single-material contracts without inventing Bi properties."""
    element = contract.get("region_material")
    if element not in ("Pb", "Bi"):
        raise ValueError("region_material must be Pb or Bi (no alloys/other materials)")
    if element == "Pb":
        expected = {"region_id": 1000, "dumpall_material_id": 2, "egs5_medium": 2}
        for key, value in expected.items():
            if contract.get(key) != value:
                raise ValueError(f"validated Pb {key} must be {value}")
        density = contract.get("cell_density_g_cm3", 9.78)
        if type(density) not in (int, float) or float(density) != 9.78:
            raise ValueError("validated Pb cell density must be 9.78 g/cm3")
        return {**expected, "region_material": "Pb", "cell_density_g_cm3": 9.78,
                "validation_status": "PB_SAME_RUN_VALIDATED"}
    for key in ("region_id", "dumpall_material_id", "egs5_medium"):
        value = contract.get(key)
        if type(value) is not int or value <= 0:
            raise ValueError(f"Bi {key} must be an explicitly configured positive integer")
    density = contract.get("cell_density_g_cm3")
    if type(density) not in (int, float) or not (float("inf") > float(density) > 0.0):
        raise ValueError("Bi cell_density_g_cm3 is required; never infer a Pb or tabulated density")
    hash_value = contract.get("stopping_table_sha256")
    if not isinstance(hash_value, str) or not re.fullmatch(r"[0-9a-fA-F]{64}", hash_value):
        raise ValueError("Bi requires the exact private HATCH stopping_table_sha256")
    return {
        "region_material": "Bi", "region_id": contract["region_id"],
        "dumpall_material_id": contract["dumpall_material_id"],
        "egs5_medium": contract["egs5_medium"],
        "cell_density_g_cm3": float(density),
        "stopping_table_sha256": hash_value.lower(),
        "validation_status": "BI_PENDING_PHITS_ORACLE",
    }


def validate_material_stopping_table(contract: dict, table_path: Path) -> str:
    """Verify bytes against the explicitly provided HATCH hash, if available.

    This only establishes file identity and material provenance; it is not
    a PHITS scoring-oracle validation of Bi.
    """
    model = material_settings(contract)
    if not table_path.is_file():
        raise FileNotFoundError(f"EGS5 HATCH table does not exist: {table_path}")
    digest = sha256(table_path)
    expected = model.get("stopping_table_sha256")
    if expected is not None and digest != expected:
        raise ValueError(f"{model['region_material']} HATCH stopping-table SHA-256 mismatch")
    return digest



def resolved_physics(settings: dict, source_input: Path) -> dict:
    """Return canonical effective model settings, excluding paths and byte hashes."""
    text = source_input.read_text(encoding="utf-8-sig", errors="replace")
    sections: dict[str, list[str]] = {}
    active = ""
    for raw in text.splitlines():
        line = raw.split("#", 1)[0].strip()
        if not line:
            continue
        if line.startswith("[") and line.endswith("]"):
            active = re.sub(r"\s+", "", line[1:-1]).lower()
            sections.setdefault(active, [])
        elif active:
            sections[active].append(line)
    return {
        "phits_version": str(settings["phits_version"]),
        "region_id": int(settings["region_id"]),
        "region_material": str(settings["region_material"]).upper(),
        "dumpall_material_id": int(settings["dumpall_material_id"]),
        "egs5_medium": int(settings["egs5_medium"]),
        "electron_cutoff_MeV": float(settings["electron_cutoff_MeV"]),
        "positron_cutoff_MeV": float(settings["positron_cutoff_MeV"]),
        "photon_cutoff_MeV": float(settings.get("photon_cutoff_MeV", 0.001)),
        "emin_sources": settings.get("emin_sources", {}),
        "esmin_MeV": float(settings["esmin_MeV"]),
        "dead_secondary_score_threshold_MeV": min(float(settings["esmin_MeV"]), 1.0e-6),
        "dead_secondary_score_rule": "PHITS-3.341-tdepstregEVENT-QCLUSTS7-gt-eletb1-v1",
        "weight_mode": str(settings["weight_mode"]).lower(),
        "energy_unit": str(settings["energy_unit"]),
        "deposit_mode": str(settings["deposit_mode"]),
        "algorithm_version": "pb-egs5-exp-v3-score-gate",
        **({"material_validation_status": "BI_PENDING_PHITS_ORACLE",
            "cell_density_g_cm3": material_settings(settings)["cell_density_g_cm3"]}
           if settings["region_material"] == "Bi" else {}),
        **source_definition(sections.get("source", [])),
    }


def physics_config(source_input: Path, config_path: Path | None = None):
    source_text = source_input.read_text(encoding="utf-8-sig", errors="replace")
    # The seed establishes run reproducibility, not local deposit physics.
    # Missing seeds must never be invented or represented as a known seed.
    seed_match = re.search(r"(?im)^\s*bitrseed\s*=\s*([01]{64})\s*(?:#.*)?$", source_text)
    sections: dict[str, list[str]] = {}
    active = ""
    for raw in source_text.splitlines():
        line = raw.split("#", 1)[0].strip()
        if not line:
            continue
        if line.startswith("[") and line.endswith("]"):
            active = re.sub(r"\s+", "", line[1:-1]).lower()
            sections.setdefault(active, [])
        elif active:
            sections[active].append(line)
    source = source_definition(sections.get("source", []))
    cutoffs = resolved_egs5_cutoffs(source_input, "3.341")
    return {
        "phits_version": "3.341", "region_id": 1000, "region_material": "PB",
        "dumpall_material_id": 2, "egs5_medium": 2,
        "electron_cutoff_MeV": cutoffs["electron_cutoff_MeV"],
        "positron_cutoff_MeV": cutoffs["positron_cutoff_MeV"],
        "photon_cutoff_MeV": cutoffs["photon_cutoff_MeV"],
        "emin_sources": cutoffs["emin_sources"],
        "esmin_MeV": resolved_esmin_MeV(source_input, "3.341"),
        "dead_secondary_score_threshold_MeV": min(resolved_esmin_MeV(source_input, "3.341"), 1.0e-6),
        "weights": "unit-weight validated slice only", "diagnostic_mode": True,
        **source,
        "unsupported_processes": "stored UNRESOLVED; never old-formula fallback",
        "source_input_sha256": sha256(source_input),
        "source_seed_bits": seed_match.group(1) if seed_match else None,
        "position_semantics": "raw transport chords and raw child point; scorer representative coordinate kept separate",
    }


def _input_section_present(path: Path, section: str) -> bool:
    """Report a PHITS section without making its settings part of the solver."""
    name = re.sub(r"\s+", "", section).lower()
    for raw in path.read_text(encoding="utf-8-sig", errors="replace").splitlines():
        candidate = raw.split("#", 1)[0].strip()
        if candidate.startswith("[") and candidate.endswith("]"):
            if re.sub(r"\s+", "", candidate[1:-1]).lower() == name:
                return True
    return False


def _input_values(path: Path, section: str) -> dict[str, str]:
    """Read the last explicit key/value occurrence in one PHITS section."""
    active = ""
    found: dict[str, str] = {}
    for raw in path.read_text(encoding="utf-8-sig", errors="replace").splitlines():
        line = raw.split("#", 1)[0].strip()
        if not line:
            continue
        if line.startswith("[") and line.endswith("]"):
            active = re.sub(r"\s+", "", line[1:-1]).lower()
        elif active == section and "=" in line:
            key, value = line.split("=", 1)
            found[re.sub(r"\s+", "", key).lower()] = value.strip()
    return found


def resolve_run_pair(dumpall: Path) -> tuple[Path, Path, str]:
    """Find exactly one sibling PHITS input/log pair for this Dumpall.

    Explicit file(15) is preferred; if omitted, use its effective value echoed
    by PHITS (or the documented dumpall.dat default). An omitted file(6) uses
    the documented phits.out default. Relative paths use the co-located-run
    convention. The referenced log must identify the executable version.
    """
    dump = dumpall.resolve(strict=True)
    if not dump.is_file():
        raise FileNotFoundError(f"Dumpall is not a regular file: {dump}")
    candidates: list[tuple[Path, Path, str]] = []
    for input_path in sorted(dump.parent.glob("*.inp")):
        params = _input_values(input_path, "parameters")
        explicit_output = params.get("file(15)")
        if explicit_output and _resolved_run_path(input_path, explicit_output) != dump:
            continue
        raw_log = params.get("file(6)", "phits.out")
        log_path = _resolved_run_path(input_path, raw_log)
        if not log_path.is_file():
            continue
        log_text = log_path.read_text(encoding="utf-8-sig", errors="replace")
        versions = set(re.findall(
            r"(?im)^\s*\|?\s*Version\s*=\s*([0-9]+\.[0-9]+)\s*\|?\s*$", log_text
        ))
        if versions != {"3.341"}:
            raise ValueError(
                f"PHITS run log {log_path} must identify exactly Version = 3.341; found {sorted(versions)}"
            )
        # PHITS echoes effective defaults in its parameter report. Use that
        # run-time value when file(15) was omitted from the source input.
        logged_output = _log_parameter(log_text, "file(15)")
        raw_output = params.get("file(15)") or logged_output or "dumpall.dat"
        output_path = _resolved_run_path(input_path, raw_output)
        if output_path != dump:
            continue
        if logged_output is not None:
            logged_path = _resolved_run_path(input_path, logged_output)
            if logged_path != output_path:
                raise ValueError(
                    f"PHITS input {input_path} file(15) resolves to {output_path}, "
                    f"but run log {log_path} reports {logged_path}"
                )
        candidates.append((input_path.resolve(), log_path, "3.341"))
    if len(candidates) != 1:
        details = ", ".join(str(item[0]) for item in candidates) or "none"
        reason = "no matching input/log pair" if not candidates else "multiple matching inputs"
        raise ValueError(
            f"cannot uniquely associate Dumpall {dump}: {reason} ({details}); "
            "expected one sibling *.inp whose explicit or PHITS-default file(15) "
            "matches the Dumpall, a file(6) or default phits.out log, and PHITS Version = 3.341"
        )
    return candidates[0]


def build_automatic_physics_contract(dumpall: Path) -> tuple[dict, Path, Path]:
    """Build the pinned Pb contract from one unambiguous same-run input/log pair."""
    source_input, run_log, version = resolve_run_pair(dumpall)
    validate_source_input(source_input)
    params = _input_values(source_input, "parameters")
    # The Dumpall stream carries NOCAS per event and is parsed by the native
    # transporter. Preserve its event order; NOCAS values need to be unique,
    # but PHITS may emit them in a non-monotonic order. Neither a maxcas count
    # nor a maxbch setting is used to compute an energy deposit. Keep present
    # settings for provenance without requiring their presence.
    seed = params.get("bitrseed", "").split()[0] if params.get("bitrseed") else None
    logged_file15 = _log_parameter(run_log.read_text(encoding="utf-8-sig", errors="replace"), "file(15)")
    effective_file15 = params.get("file(15)") or logged_file15 or "dumpall.dat"
    maxcas_raw = params.get("maxcas", "").split()[:1]
    maxbch_raw = params.get("maxbch", "").split()[:1]
    try:
        maxcas_decimal = Decimal(maxcas_raw[0].replace("D", "E").replace("d", "e")) if maxcas_raw else None
        maxbch_decimal = Decimal(maxbch_raw[0].replace("D", "E").replace("d", "e")) if maxbch_raw else None
        if (maxcas_decimal is not None and
                (not maxcas_decimal.is_finite() or maxcas_decimal != maxcas_decimal.to_integral_value())):
            raise ValueError("maxcas is not an integer-valued number")
        if (maxbch_decimal is not None and
                (not maxbch_decimal.is_finite() or maxbch_decimal != maxbch_decimal.to_integral_value())):
            raise ValueError("maxbch is not an integer-valued number")
        maxcas = int(maxcas_decimal) if maxcas_decimal is not None else None
        maxbch = int(maxbch_decimal) if maxbch_decimal is not None else None
    except (InvalidOperation, ValueError) as error:
        raise ValueError("explicit maxcas/maxbch must be integer-valued PHITS numbers") from error
    if maxcas is not None and maxcas <= 0:
        raise ValueError("explicit maxcas must be positive")
    if maxbch is not None and maxbch <= 0:
        raise ValueError("explicit maxbch must be positive")
    esmin = resolved_esmin_MeV(source_input, version)
    cutoffs = resolved_egs5_cutoffs(source_input, version)
    contract = {
        "phits_version": version,
        "phits_input": str(source_input),
        "phits_input_sha256": sha256(source_input),
        "phits_run_log": str(run_log),
        "phits_run_log_sha256": sha256(run_log),
        "phits_file6_status": "explicit" if params.get("file(6)") else "phits_default",
        "phits_file15_status": ("explicit" if params.get("file(15)") else
                                 "phits_log_effective" if logged_file15 is not None else "phits_default"),
        "phits_file15_effective": str(_resolved_run_path(source_input, effective_file15)),
        "source_seed_bits": seed,
        "source_seed_status": "explicit" if seed else "not_explicit_in_input",
        "maxcas": maxcas,
        "maxbch": maxbch,
        "maxcas_status": "explicit" if maxcas is not None else "not_explicit_in_input",
        "maxbch_status": "explicit" if maxbch is not None else "not_explicit_in_input",
        "region_id": 1000,
        "region_material": "Pb",
        "dumpall_material_id": 2,
        "egs5_medium": 2,
        "electron_cutoff_MeV": cutoffs["electron_cutoff_MeV"],
        "positron_cutoff_MeV": cutoffs["positron_cutoff_MeV"],
        "photon_cutoff_MeV": cutoffs["photon_cutoff_MeV"],
        "emin_sources": cutoffs["emin_sources"],
        "esmin_MeV": esmin,
        "dead_secondary_score_threshold_MeV": min(esmin, 1.0e-6),
        "weight_mode": "unit",
        "energy_unit": "MeV",
        "deposit_mode": "t-deposit-output-deposit",
    }
    validate_source_input(source_input, contract)
    return contract, source_input, run_log


def provenance_metadata(physics_settings: dict | None, config_path: Path | None,
                        source_input: Path, dumpall: Path, stopping_table: Path,
                        converter_executable: Path | None = None,
                        reconstructor: Path | None = None,
                        track_exporter: Path | None = None) -> dict:
    """Build distinct byte-level and semantic hashes for the HDF5 v2 provenance."""
    settings = physics_settings or {
        "phits_version": "3.341", "region_id": 1000, "region_material": "Pb",
        "dumpall_material_id": 2, "egs5_medium": 2, "electron_cutoff_MeV": 0.001,
        "positron_cutoff_MeV": 0.001, "weight_mode": "unit", "energy_unit": "MeV",
        "esmin_MeV": 1.0e-3,
        "deposit_mode": "t-deposit-output-deposit",
    }
    resolved = resolved_physics(settings, source_input)
    resolved_hash = hashlib.sha256(canonical_json(resolved)).hexdigest()
    input_hash = sha256(source_input)
    physics = physics_config(source_input)
    metadata = {
        "provenance_schema_revision": 2,
        # Kept as a v2 compatibility alias, now defined identically to resolved_physics_sha256.
        "physics_config_sha256": resolved_hash,
        "resolved_physics_sha256": resolved_hash,
        "phits_input_sha256": input_hash,
        "source_input_sha256": input_hash,
        "dumpall_sha256": sha256(dumpall),
        "stopping_table_sha256": sha256(stopping_table),
        "phits_version": str(settings["phits_version"]),
        "resolved_physics_json": canonical_json(resolved).decode("utf-8"),
        "source_mode": physics["source_mode"],
        "source_energy_MeV": ",".join(format(x, ".12g") for x in physics["source_energies_MeV"]),
        "reconstruction_algorithm_version": "pb-egs5-exp-v3-score-gate",
        "validation_scope": "PHITS-3.341; Pb region 1000; unit-weight; validated EGS5 process subset",
        "tally_semantics": "Dumpall reconstructed locally; PHITS [t-deposit] is optional and not read; transport-step position remains a segment",
        "t_deposit_status": "configured_not_used" if _input_section_present(source_input, "t-deposit") else "absent",
        "source_seed_status": "explicit" if physics["source_seed_bits"] else "not_explicit_in_input",
        "driver_sha256": sha256(Path(__file__).resolve()),
    }
    if config_path is not None:
        metadata["physics_contract_sha256"] = sha256(config_path)
        contract = json.loads(config_path.read_text(encoding="utf-8"))
        phits_executable_sha = contract.get("phits_executable_sha256")
        if phits_executable_sha is not None:
            if not re.fullmatch(r"[0-9a-fA-F]{64}", str(phits_executable_sha)):
                raise ValueError("phits_executable_sha256 in physics contract must be a 64-character SHA-256 digest")
            metadata["phits_executable_sha256"] = str(phits_executable_sha).lower()
    for key, path in (("converter_executable_sha256", converter_executable),
                      ("reconstructor_executable_sha256", reconstructor),
                      ("transport_exporter_sha256", track_exporter)):
        if path is not None:
            metadata[key] = sha256(path)
    return metadata


def validate_source_input(path: Path, contract: dict | None = None) -> None:
    """Validate a one-material PHITS scoring input without asserting Bi oracle closure."""
    model = material_settings(contract or {"region_material": "Pb", "region_id": 1000,
                                                "dumpall_material_id": 2, "egs5_medium": 2})
    text = path.read_text(encoding="utf-8-sig", errors="replace")
    active = ""
    sections: dict[str, list[str]] = {}
    for raw in text.splitlines():
        line = raw.split("#", 1)[0].strip()
        if not line:
            continue
        if line.startswith("[") and line.endswith("]"):
            active = re.sub(r"\s+", "", line[1:-1]).lower()
            sections.setdefault(active, [])
        elif active:
            sections[active].append(line.lower())
    def values(section: str) -> dict[str, str]:
        found = {}
        for line in sections.get(section, []):
            if "=" not in line:
                continue
            key, value = line.split("=", 1)
            found[re.sub(r"\s+", "", key)] = value.strip()
        return found
    params = values("parameters")
    # Resolve omitted settings to PHITS defaults before enforcing the
    # converter's supported mode. A different effective cutoff is not silently
    # treated as the validated 1 keV setting.
    mode_values = {"negs": ("1", "-1"), "e-mode": ("2", "0"),
                   "icntl": ("0", "0"), "dumpall": ("-1", "0")}
    for key, (expected, default) in mode_values.items():
        raw = params.get(key)
        actual = raw.split()[0] if raw is not None else default
        source = "explicit" if raw is not None else "PHITS default"
        if _parameter_number(actual, key) != float(expected):
            raise ValueError(
                f"effective [{key}]={actual} ({source}); this conversion requires {expected}"
            )
    cutoffs = resolved_egs5_cutoffs(path, "3.341")
    cutoff_keys = {12: "electron_cutoff_MeV", 13: "positron_cutoff_MeV", 14: "photon_cutoff_MeV"}
    for index, expected in VALIDATED_EGS5_CUTOFFS_MEV.items():
        key = f"emin({index})"
        actual = cutoffs[cutoff_keys[index]]
        source = cutoffs["emin_sources"][key]
        if actual != expected:
            raise ValueError(
                f"effective [{key}]={actual:g} MeV ({source}); the validated Pb "
                f"reconstruction requires {expected:g} MeV"
            )
    seed = params.get("bitrseed", "").split()[0] if params.get("bitrseed") else None
    if seed is not None and not re.fullmatch(r"[01]{64}", seed):
        raise ValueError("explicit bitrseed must contain exactly 64 binary digits")
    material_lines = sections.get("material", [])
    material_id = model["dumpall_material_id"]
    element = model["region_material"].lower()
    expected_material = rf"mat\s*\[\s*{material_id}\s*\]\s+{element}\s+1(?:\s|$)"
    if not any(re.match(expected_material, line) for line in material_lines):
        raise ValueError(f"source input must define mat[{material_id}] as pure {model['region_material']}")
    cell_lines = sections.get("cell", [])
    region = model["region_id"]
    density = model["cell_density_g_cm3"]
    # Match only the region's single-cell definition, never an unrelated material.
    cell_matches = [
        re.match(rf"^{region}\s+(\d+)\s+(-?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?)(?:\s|$)", line)
        for line in cell_lines if re.match(rf"^{region}(?:\s|$)", line)
    ]
    if len(cell_matches) != 1 or cell_matches[0] is None:
        raise ValueError(f"source input must define exactly one simple cell for region {region}")
    found_material = int(cell_matches[0].group(1))
    found_density = float(cell_matches[0].group(2))
    if found_material != material_id or found_density != -density:
        raise ValueError(f"region {region} requires material {material_id} and density -{density:g}")
    source_definition(sections.get("source", []))
    esmin = resolved_esmin_MeV(path, "3.341")
    # Force evaluation now; the exact score gate is captured in physics config and provenance.
    if not (esmin > 0.0 and esmin < float("inf")):
        raise ValueError("PHITS 3.341 esmin did not resolve to a finite positive value")
    # [t-deposit] affects independent PHITS tally output, NOT the Dumpall
    # records or our reconstruction. No tally section/binning is required.
    # Never treat its absence (or a different histogram) as a failure of the
    # material/cutoff/Dumpall transport contract.


def load_physics_config(path: Path, *, dumpall: Path | None = None,
                        allow_pending_bi: bool = False) -> tuple[dict, Path]:
    data = json.loads(path.read_text(encoding="utf-8"))
    model = material_settings(data)
    expected = {
        "phits_version": "3.341",
        "electron_cutoff_MeV": 0.001, "positron_cutoff_MeV": 0.001,
        "weight_mode": "unit", "energy_unit": "MeV", "deposit_mode": "t-deposit-output-deposit",
    }
    for key, value in expected.items():
        if data.get(key) != value:
            raise ValueError(f"physics config {key} must be {value!r}; found {data.get(key)!r}")
    source_value = data.get("phits_input")
    if not source_value:
        raise ValueError("physics config must name the PHITS input in 'phits_input'")
    source_input = Path(source_value)
    if not source_input.is_absolute():
        source_input = (path.parent / source_input).resolve()
    if not source_input.is_file():
        raise FileNotFoundError(f"PHITS input from physics config not found: {source_input}")
    if dumpall is not None:
        paired_input, _run_log, run_version = resolve_run_pair(dumpall)
        if paired_input != source_input.resolve():
            raise ValueError(
                f"physics config PHITS input {source_input.resolve()} is not the unique input "
                f"associated with Dumpall {dumpall.resolve()} ({paired_input})"
            )
        if run_version != data["phits_version"]:
            raise ValueError(
                f"physics config PHITS version {data['phits_version']!r} does not match run log {run_version!r}"
            )
    validate_source_input(source_input, data)
    if model["region_material"] == "Bi" and not allow_pending_bi:
        raise ValueError("Bi physics is not yet PHITS-oracle validated; HDF5 publication is disabled")
    resolved_esmin = resolved_esmin_MeV(source_input, str(data["phits_version"]))
    if float(data.get("esmin_MeV", float("nan"))) != resolved_esmin:
        raise ValueError(f"physics config esmin_MeV {data.get('esmin_MeV')!r} does not match resolved PHITS value {resolved_esmin}")
    expected_threshold = min(resolved_esmin, 1.0e-6)
    if float(data.get("dead_secondary_score_threshold_MeV", float("nan"))) != expected_threshold:
        raise ValueError("physics config dead_secondary_score_threshold_MeV does not equal min(esmin,1e-6)")
    return data, source_input


def groups(path: Path):
    stream = path.open("r", encoding="utf-8-sig", newline="")
    reader = csv.DictReader(stream)
    if not reader.fieldnames or "event_id" not in reader.fieldnames:
        stream.close()
        raise ValueError(f"{path} lacks event_id")
    grouped = itertools.groupby(reader, key=lambda row: int(row["event_id"]))
    try:
        for event_id, rows in grouped:
            yield event_id, list(rows)
    finally:
        stream.close()


def take(group_iter, current, event_id):
    if current is None:
        current = next(group_iter, None)
    if current is None:
        return [], None
    if current[0] != event_id:
        raise ValueError(
            f"input event stream is not aligned with raw event order: "
            f"expected event {event_id}, found {current[0]}"
        )
    return current[1], next(group_iter, None)


def run(args, writer: EventV2Writer, provenance: dict) -> dict:
    for path in (args.dumpall, args.stopping_table, args.reconstructor, args.track_exporter):
        if not path.is_file():
            raise FileNotFoundError(path)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.work_dir.mkdir(parents=True, exist_ok=True)
    ledger_path = args.work_dir / "physical_ledger.csv"
    not_scored_path = args.work_dir / "not_scored_candidates.csv"
    event_path = args.work_dir / "raw_events.csv"
    track_path = args.work_dir / "raw_tracks.csv"
    threshold = float(args.dead_secondary_score_threshold_MeV)
    ledger_run = subprocess.run([str(args.reconstructor), str(args.dumpall),
                                 str(args.stopping_table), str(ledger_path),
                                 str(not_scored_path), format(threshold, ".17g")],
                                check=False, text=True, capture_output=True)
    if ledger_run.returncode not in (0, 1) or not ledger_path.is_file():
        raise RuntimeError(f"physical reconstruction failed: rc={ledger_run.returncode}\n{ledger_run.stderr}")
    parse_run = subprocess.run([str(args.track_exporter), str(args.dumpall),
                                 str(event_path), str(track_path)],
                               check=False, text=True, capture_output=True)
    if parse_run.returncode != 0:
        raise RuntimeError(f"raw transport export failed: rc={parse_run.returncode}\n{parse_run.stderr}")
    with not_scored_path.open("r", encoding="utf-8-sig", newline="") as handle:
        not_scored_rows = list(csv.DictReader(handle))
    if any(row.get("classification") != "NOT_SCORED_BY_PHITS" for row in not_scored_rows):
        raise ValueError("non-scored candidate audit contains an unexpected classification")

    track_groups = iter(groups(track_path))
    ledger_groups = iter(groups(ledger_path))
    next_tracks = next(track_groups, None)
    next_ledger = next(ledger_groups, None)
    event_count = track_count = particle_count = deposit_count = unresolved_count = 0
    negative_count = duplicate_ids = 0
    seen_event_ids = set()
    with event_path.open("r", encoding="utf-8-sig", newline="") as event_stream:
        for row in csv.DictReader(event_stream):
            event_id = int(row["event_id"])
            if event_id <= 0 or event_id in seen_event_ids:
                raise ValueError("raw event IDs must be positive and unique")
            seen_event_ids.add(event_id)
            raw_tracks, next_tracks = take(track_groups, next_tracks, event_id)
            raw_deposits, next_ledger = take(ledger_groups, next_ledger, event_id)

            particles_by_id = {}
            by_source_record = {}
            for track in raw_tracks:
                pid = int(track["particle_id"])
                ityp = int(track["ityp"])
                particle = particles_by_id.setdefault(pid, {"particle_id": pid,
                    "parent_particle_id": None, "ityp": ityp, "tracks": []})
                if particle["ityp"] != ityp:
                    raise ValueError(f"particle {pid} changes ITYP in event {event_id}")
                tid = track["source_record_id"]
                parsed_track = {"source_record_id": tid, "ncol": int(track["ncol"]),
                    "xyz_cm": (float(track["x_cm"]), float(track["y_cm"]), float(track["z_cm"])),
                    "energy_MeV": float(track["energy_MeV"]), "time": float(track["time"]),
                    "region_id": int(track["region_id"])}
                particle["tracks"].append(parsed_track)
                by_source_record[tid] = parsed_track
            particles = list(particles_by_id.values())
            deposits = []
            event_unresolved = 0
            for raw in raw_deposits:
                status = raw["status"]
                energy = float(raw["E_deposit_MeV"]) if raw["E_deposit_MeV"] else None
                negative_count += int(energy is not None and energy < 0)
                kind = "UNRESOLVED" if status != "reconstructed" else "SEGMENT" if raw["position_status"].startswith("segment") else "POINT"
                pos = None
                if kind == "SEGMENT":
                    start = tuple(float(raw[k]) for k in ("x", "y", "z"))
                    end = tuple(float(raw[k]) for k in ("x_end", "y_end", "z_end"))
                    score = None
                    semantics = "EGS5 continuous loss interval from raw Dumpall endpoints; score representative point unavailable"
                elif kind == "POINT" and raw["position_status"] == "secondary_point":
                    pos = tuple(float(raw[k]) for k in ("x", "y", "z"))
                    start = end = None
                    score = pos
                    semantics = "raw child position; local point validated for dead charged e-/e+ secondary branches"
                elif kind == "POINT" and energy == 0.0:
                    record = by_source_record.get(raw["source_line"])
                    if record is None:
                        raise ValueError(f"zero addition {raw['addition_id']} has no raw transport point")
                    pos = record["xyz_cm"]
                    start = end = None
                    score = pos
                    semantics = "raw transport point retained for a zero-energy tally addition; no heat contribution"
                else:
                    start = end = score = None
                    semantics = "unresolved position"
                deposits.append({
                    "addition_id": raw["addition_id"], "particle_id": int(raw["particle_id"]),
                    "process": raw["process"], "region_id": int(raw["region_id"]),
                    "energy_MeV": energy, "kind": kind, "start_xyz_cm": start,
                    "end_xyz_cm": end, "score_xyz_cm": score,
                    "source_record_id": raw["source_line"], "position_semantics": semantics,
                })
                deposit_count += 1
                unresolved_count += int(kind == "UNRESOLVED")
                event_unresolved += int(kind == "UNRESOLVED")
            track_count += len(raw_tracks)
            particle_count += len(particles)
            complete = event_unresolved == 0
            writer.append_event(event_id, particles, deposits, complete=complete,
                                all_absorption="undetermined")
            event_count += 1
    if next_tracks is not None or next_ledger is not None:
        raise ValueError("track or deposit input contains events absent from the raw event stream")
    if negative_count:
        raise ValueError(f"negative physical energies found: {negative_count}")
    if unresolved_count:
        raise ValueError(f"unsupported scoring branches remain ({unresolved_count}); refusing to publish complete HDF5 v2")
    return {
        "input_dumpall": str(args.dumpall), "input_bytes": args.dumpall.stat().st_size,
        "dumpall_sha256": sha256(args.dumpall), "stopping_table_sha256": sha256(args.stopping_table),
        **provenance,
        "physical_ledger_exit_code": ledger_run.returncode,
        "physical_ledger_stdout": ledger_run.stdout.strip(), "transport_export_stdout": parse_run.stdout.strip(),
        "not_scored_by_phits_candidates": len(not_scored_rows),
        "not_scored_candidate_audit": str(not_scored_path),
        "events": event_count, "particles": particle_count, "tracks": track_count,
        "deposits": deposit_count, "unresolved": unresolved_count, "negative": negative_count,
        "output": str(args.output),
    }


# Pinned to the exact, same-run-PHITS-validated historical HATCH CSV.
# This data file is NOT automatically published here: verify distribution
# permission before adding it to a public Git repository.
BUNDLED_PB_TABLE_NAME = "pb_egs5_phits3341.csv"
BUNDLED_PB_TABLE_SHA256 = "20fd86eb0b57483bd19f20a8388d2c1ecb0431a02fa99f0933e2c3c4748514e9"


def resolve_stopping_table(explicit: Path | None) -> Path:
    """Select an explicitly provided table or an exact, checksum-pinned bundled Pb table.

    No unverified physics table is silently substituted, and no synthetic
    fixture is ever a valid default.
    """
    if explicit is not None:
        candidates = [explicit]
    else:
        # Prefer the immutable, SHA-pinned distribution resource over any old
        # per-user environment variable. The environment path is compatibility
        # fallback only when no bundled table is installed.
        home = Path(__file__).resolve().parent
        candidates = [
            home / BUNDLED_PB_TABLE_NAME,  # Python package install beside driver
            home.parent / "data" / BUNDLED_PB_TABLE_NAME,  # repository source tree
            home.parent.parent / "tes_cpp" / "bin" / BUNDLED_PB_TABLE_NAME,  # scikit-build/CMake install prefix
            # CMake/scikit-build installs native assets in the Python package
            # directory while placing this driver under share/.
            home.parent.parent / "tes_cpp_dump2event" / "bin" / BUNDLED_PB_TABLE_NAME,
            home.parent.parent / "bin" / BUNDLED_PB_TABLE_NAME,  # CMake prefix/bin
        ]
        configured = os.environ.get("TES_CPP_PB_HATCH_TABLE")
        if configured:
            candidates.append(Path(configured).expanduser())
    for candidate in candidates:
        if candidate.is_file():
            actual = sha256(candidate)
            if actual != BUNDLED_PB_TABLE_SHA256:
                raise ValueError(
                    f"Pb HATCH coefficient checksum mismatch: {candidate} "
                    f"(expected {BUNDLED_PB_TABLE_SHA256}, got {actual})"
                )
            return candidate
    if explicit is not None:
        raise FileNotFoundError(f"explicit stopping table was not found: {explicit}")
    raise FileNotFoundError(
        "verified Pb HATCH coefficient table is not bundled. "
        "Install the byte-verified table as tes_cpp/dump2event/data/"
        + BUNDLED_PB_TABLE_NAME +
        "; set TES_CPP_PB_HATCH_TABLE for a private installation, or specify "
        "--stopping-table with the exact SHA-256-pinned CSV."
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dumpall", type=Path, required=True)
    parser.add_argument("--stopping-table", type=Path,
                        help="optional override; omitted means the checksum-pinned bundled Pb HATCH CSV")
    parser.add_argument("--source-input", type=Path,
                        help="legacy: PHITS input whose settings are validated")
    parser.add_argument("--physics-config", type=Path,
                        help="explicit JSON contract for the experimental Pb EGS5 mode")
    parser.add_argument("--reconstructor", type=Path, required=True)
    parser.add_argument("--track-exporter", type=Path, required=True)
    parser.add_argument("--converter-executable", type=Path,
                        help="the dump2event executable that launched this driver")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--work-dir", type=Path,
                        help="optional scratch directory (otherwise a private temporary directory is used)")
    args = parser.parse_args()
    try:
        if args.output.exists():
            raise FileExistsError(f"refusing to overwrite existing output: {args.output}")
        args.stopping_table = resolve_stopping_table(args.stopping_table)
        if args.physics_config and args.source_input:
            raise ValueError("provide only one of --physics-config and --source-input")
        automatic_contract = None
        if args.physics_config:
            physics_data, args.source_input = load_physics_config(
                args.physics_config, dumpall=args.dumpall)
        elif args.source_input:
            if not args.source_input.is_file():
                raise FileNotFoundError(args.source_input)
            # Follow the same canonical validated contract as automatic mode.
            # Previously this branch left physics_data=None and crashed when
            # building the HDF5 material provenance.
            paired_contract, paired_input, _run_log = build_automatic_physics_contract(args.dumpall)
            if paired_input != args.source_input.resolve():
                raise ValueError("--source-input does not match the unique PHITS input associated with Dumpall")
            args.source_input = paired_input
            automatic_contract = paired_contract
            physics_data = paired_contract
        else:
            automatic_contract, args.source_input, run_log = build_automatic_physics_contract(args.dumpall)
            physics_data = automatic_contract
        # The writer is held open only for this conversion. Output publication is
        # atomic; input data and existing destinations are never modified.
        physics = physics_config(args.source_input)
        args.dead_secondary_score_threshold_MeV = float(
            physics_data["dead_secondary_score_threshold_MeV"] if physics_data else
            physics["dead_secondary_score_threshold_MeV"])
        provenance = provenance_metadata(physics_data, args.physics_config,
                                         args.source_input, args.dumpall, args.stopping_table,
                                         args.converter_executable, args.reconstructor, args.track_exporter)
        if automatic_contract is not None:
            automatic_json = canonical_json(automatic_contract)
            automatic_hash = hashlib.sha256(automatic_json).hexdigest()
            provenance.update({
                "physics_config_origin": "auto-associated-same-folder-phits-input",
                "physics_contract_json": automatic_json.decode("utf-8"),
                "physics_contract_sha256": automatic_hash,
                "generated_physics_config_json": automatic_json.decode("utf-8"),
                "generated_physics_config_sha256": automatic_hash,
                "phits_run_log_sha256": automatic_contract["phits_run_log_sha256"],
            })
        config = {**provenance, "material": physics_data["region_material"],
                  "medium": str(physics_data["egs5_medium"]), "cutoff_MeV": "0.001",
                  "algorithm_version": "pb-egs5-exp-v3-score-gate",
                  "position_semantics": physics["position_semantics"],
                  "unsupported_processes": physics["unsupported_processes"],
                  "weights": physics["weights"], "diagnostic_mode": True,
                  "tally_semantics": "Dumpall reconstructed locally; [t-deposit] optional and not used; transport-step position remains a segment",
                  "complete_semantics": "per-event scoring complete only when no in-scope scored candidate is UNRESOLVED; proven not-scored candidates are excluded from deposits and preserved in audit",
                  "dead_secondary_score_rule": "PHITS-3.341-tdepstregEVENT-QCLUSTS7-gt-eletb1-v1",
                  "dead_secondary_score_threshold_MeV": format(float(physics_data["dead_secondary_score_threshold_MeV"] if physics_data else physics["dead_secondary_score_threshold_MeV"]), ".17g"),
                  "esmin_MeV": format(float(physics_data["esmin_MeV"] if physics_data else physics["esmin_MeV"]), ".17g"),
                  "dumpall_structure_semantics": "conversion requires terminal NCOL=2",
                  "transport_physics_completeness": "undetermined",
                  "energy_balance_status": "not assessed",
                  "all_absorption": "undetermined",
                  "region_id": physics["region_id"], "dumpall_material_id": physics["dumpall_material_id"],
                  "region_material": physics["region_material"], "egs5_medium": physics["egs5_medium"]}
        started = time.perf_counter()
        # load_physics_config refuses Bi until an independent same-run PHITS
        # scoring oracle validates that material and source configuration.
        with tempfile.TemporaryDirectory(prefix="dump2event-egs5-", dir=args.output.parent) as scratch:
            if args.work_dir:
                args.work_dir.mkdir(parents=True, exist_ok=True)
            else:
                args.work_dir = Path(scratch)
            with EventV2Writer(args.output, config) as writer:
                result = run(args, writer, provenance)
        result["output_bytes"] = args.output.stat().st_size
        result["elapsed_seconds"] = time.perf_counter() - started
        print(json.dumps(result, indent=2))
        return 0
    except Exception as error:
        print(f"INVALID: {error}")
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
