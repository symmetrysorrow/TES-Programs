"""Vendored 2026-07-14 from Thermal-and-Electoric-Sim core/geometry/body_semantics.py
(import paths adjusted only). See vendored/__init__.py for policy."""
from __future__ import annotations


def body_name_of(item: object) -> str:
    return (
        str(getattr(item, "body_name", "") or "").strip()
        or str(getattr(item, "group_name", "") or "").strip()
        or str(getattr(item, "name", "") or "").strip()
    )


def explicit_body_name_of(item: object) -> str:
    return str(getattr(item, "body_name", "") or "").strip()


def body_mode_of(item: object) -> str:
    return (
        str(getattr(item, "body_mode", "") or "").strip().lower()
        or str(getattr(item, "group_mode", "") or "").strip().lower()
        or "add"
    )


def uses_body_boolean(item: object) -> bool:
    explicit_body_name = explicit_body_name_of(item)
    return body_mode_of(item) != "add" or (explicit_body_name and explicit_body_name != str(getattr(item, "name", "") or "").strip())
