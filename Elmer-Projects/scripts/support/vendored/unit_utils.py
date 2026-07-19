"""Physical unit/dimension tables (vendored subset of core/project/unit_utils.py).

Vendored 2026-07-14 from Thermal-and-Electoric-Sim; trimmed to what
dimensioned_expression.py needs (display helpers and GUI-side utilities
removed). Do not extend here — sync from upstream instead.
"""
from __future__ import annotations

from typing import Any

_UNIT_CHAR_ALIASES: dict[str, str] = {"µ": "u", "μ": "u", "Ω": "ohm"}
_UNIT_NAME_ALIASES: dict[str, str] = {
    "w/k": "wperkelvin",
    "mw/k": "mwperkelvin",
    "uw/k": "uwperkelvin",
    "nw/k": "nwperkelvin",
    "w/(m*k)": "wpermeterkelvin",
    "mw/(m*k)": "mwpermeterkelvin",
    "uw/(m*k)": "uwpermeterkelvin",
    "nw/(m*k)": "nwpermeterkelvin",
    "w/(m^2*k)": "wpermetersquaredkelvin",
    "mw/(m^2*k)": "mwpermetersquaredkelvin",
    "uw/(m^2*k)": "uwpermetersquaredkelvin",
    "nw/(m^2*k)": "nwpermetersquaredkelvin",
    "w/(m^2)": "wpermetersquared",
    "mw/(m^2)": "mwpermetersquared",
    "uw/(m^2)": "uwpermetersquared",
    "nw/(m^2)": "nwpermetersquared",
    "kg/(m^3)": "kgpermetercubed",
    "g/(m^3)": "gpermetercubed",
    "mg/(m^3)": "mgpermetercubed",
    "ug/(m^3)": "ugpermetercubed",
    "j/(kg*k)": "jperkgkelvin",
    "kj/(kg*k)": "kjperkgkelvin",
    "mj/(kg*k)": "mjperkgkelvin",
    "uj/(kg*k)": "ujperkgkelvin",
    "(m^2*k)/w": "squaremeterkelvinperwatt",
    "(mm^2*k)/w": "millisquaremeterkelvinperwatt",
    "(um^2*k)/w": "microsquaremeterkelvinperwatt",
    "1/k": "perkelvin",
}


def normalize_unit(unit: str) -> str:
    """Normalize micro-sign variants and lower-case a unit string."""
    s = str(unit or "").strip()
    for src, dst in _UNIT_CHAR_ALIASES.items():
        s = s.replace(src, dst)
    normalized = s.lower().replace(" ", "")
    return _UNIT_NAME_ALIASES.get(normalized, normalized)


_BASE_DIMENSIONS: tuple[str, ...] = ("mass", "length", "time", "current", "temperature")


def _dims(**powers: int) -> tuple[tuple[str, int], ...]:
    return _normalize_signature(
        (name, int(powers[name])) for name in _BASE_DIMENSIONS if int(powers.get(name, 0)) != 0
    )


def _normalize_signature(signature: Any) -> tuple[tuple[str, int], ...]:
    merged: dict[str, int] = {}
    for dim, power in signature:
        dim_name = str(dim).strip().lower()
        power_i = int(power)
        if power_i == 0:
            continue
        merged[dim_name] = merged.get(dim_name, 0) + power_i
    return tuple(sorted((dim, power) for dim, power in merged.items() if power != 0))


_DIMENSION_SIGNATURES: dict[str, tuple[tuple[str, int], ...]] = {
    "dimensionless": (),
    "length": _dims(length=1),
    "time": _dims(time=1),
    "energy": _dims(mass=1, length=2, time=-2),
    "temperature": _dims(temperature=1),
    "resistance": _dims(mass=1, length=2, time=-3, current=-2),
    "inductance": _dims(mass=1, length=2, time=-2, current=-2),
    "voltage": _dims(mass=1, length=2, time=-3, current=-1),
    "current": _dims(current=1),
    "pressure": _dims(mass=1, length=-1, time=-2),
    "conductivity": _dims(mass=1, length=1, time=-3, temperature=-1),
    "thermal_conductance": _dims(mass=1, length=2, time=-3, temperature=-1),
    "heat_transfer_coefficient": _dims(mass=1, time=-3, temperature=-1),
    "heat_flux": _dims(mass=1, time=-3),
    "density": _dims(mass=1, length=-3),
    "specific_heat": _dims(length=2, time=-2, temperature=-1),
    "area_thermal_resistance": _dims(mass=-1, time=3, temperature=1),
    "per_kelvin": _dims(temperature=-1),
}

# NOTE (upstream domain convention): the energy base unit is keV, not J.
_PHYSICAL_DIMENSIONS: dict[str, dict[str, float]] = {
    "length":      {"km": 1e3,  "m": 1.0,  "mm": 1e-3,  "um": 1e-6,  "nm": 1e-9,  "pm": 1e-12},
    "time":        {"s": 1.0,   "ms": 1e-3, "us": 1e-6,  "ns": 1e-9,  "ps": 1e-12},
    "energy":      {"mev": 1e3, "kev": 1.0, "ev": 1e-3, "j": 6.241509074460763e15},
    "temperature": {"k": 1.0,  "mk": 1e-3},
    "resistance":  {"kohm": 1e3, "ohm": 1.0, "mohm": 1e-3},
    "inductance":  {"h": 1.0,  "mh": 1e-3, "uh": 1e-6, "nh": 1e-9, "ph": 1e-12},
    "voltage":     {"v": 1.0,  "mv": 1e-3, "uv": 1e-6},
    "current":     {"a": 1.0,  "ma": 1e-3, "ua": 1e-6, "na": 1e-9},
    "pressure":    {"pa": 1.0, "kpa": 1e3, "mpa": 1e6, "gpa": 1e9},
    "conductivity": {
        "wpermeterkelvin": 1.0,
        "mwpermeterkelvin": 1e-3,
        "uwpermeterkelvin": 1e-6,
        "nwpermeterkelvin": 1e-9,
    },
    "thermal_conductance": {
        "wperkelvin": 1.0,
        "mwperkelvin": 1e-3,
        "uwperkelvin": 1e-6,
        "nwperkelvin": 1e-9,
    },
    "heat_transfer_coefficient": {
        "wpermetersquaredkelvin": 1.0,
        "mwpermetersquaredkelvin": 1e-3,
        "uwpermetersquaredkelvin": 1e-6,
        "nwpermetersquaredkelvin": 1e-9,
    },
    "heat_flux": {
        "wpermetersquared": 1.0,
        "mwpermetersquared": 1e-3,
        "uwpermetersquared": 1e-6,
        "nwpermetersquared": 1e-9,
    },
    "density": {
        "kgpermetercubed": 1.0,
        "gpermetercubed": 1e-3,
        "mgpermetercubed": 1e-6,
        "ugpermetercubed": 1e-9,
    },
    "specific_heat": {
        "jperkgkelvin": 1.0,
        "kjperkgkelvin": 1e3,
        "mjperkgkelvin": 1e-3,
        "ujperkgkelvin": 1e-6,
    },
    "area_thermal_resistance": {
        "squaremeterkelvinperwatt": 1.0,
        "millisquaremeterkelvinperwatt": 1e-3,
        "microsquaremeterkelvinperwatt": 1e-6,
    },
    "per_kelvin": {
        "perkelvin": 1.0,
        "mperkelvin": 1e-3,
        "uperkelvin": 1e-6,
        "nperkelvin": 1e-9,
    },
    "dimensionless": {},
}

_UNIT_TO_DIMENSION: dict[str, str] = {
    unit: dim
    for dim, units in _PHYSICAL_DIMENSIONS.items()
    for unit in units
}

_SIGNATURE_TO_DIMENSION: dict[tuple[tuple[str, int], ...], str] = {
    signature: name
    for name, signature in _DIMENSION_SIGNATURES.items()
    if signature
}


def dimension_signature_for_name(name: str) -> tuple[tuple[str, int], ...] | None:
    normalized = normalize_unit(name)
    if normalized in _DIMENSION_SIGNATURES:
        return _DIMENSION_SIGNATURES[normalized]
    if normalized in _BASE_DIMENSIONS:
        return ((normalized, 1),)
    return None


def dimension_name_from_signature(signature: tuple[tuple[str, int], ...]) -> str | None:
    normalized = _normalize_signature(signature)
    if not normalized:
        return None
    return _SIGNATURE_TO_DIMENSION.get(normalized)


def unit_factor_and_dimension(unit: str) -> tuple[tuple[tuple[str, int], ...], float] | None:
    norm = normalize_unit(unit)
    dim = _UNIT_TO_DIMENSION.get(norm)
    if dim is None:
        return None
    factor = _PHYSICAL_DIMENSIONS.get(dim, {}).get(norm)
    if factor is None:
        return None
    signature = _DIMENSION_SIGNATURES.get(dim)
    if signature is None:
        return None
    return signature, float(factor)


def units_for_dimension(dimension: str) -> dict[str, float]:
    """Return a ``{normalised_unit: SI_factor}`` dict for *dimension*."""
    return dict(_PHYSICAL_DIMENSIONS.get(dimension, {}))
