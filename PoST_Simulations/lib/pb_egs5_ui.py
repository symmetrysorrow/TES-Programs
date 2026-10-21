"""Interactive, opt-in PHITS Pb EGS5 -> HDF5 v2 conversion.

This UI intentionally does not call Pulse_Ms.  The v2 deposit schema has
SEGMENT records that require a separately validated TES allocation model.
"""

from __future__ import annotations

import json
from pathlib import Path


def _ask_path(questionary, message: str, *, default: str = "") -> Path | None:
    answer = questionary.path(message, default=default).ask()
    if answer is None:
        return None
    if not str(answer).strip():
        raise ValueError("ファイルパスが空です")
    return Path(str(answer).strip()).expanduser()


def _checked_input(path: Path, name: str) -> Path:
    if not path.is_file():
        raise FileNotFoundError(f"{name}が見つかりません: {path}")
    return path.resolve()


def _read_contract_input(contract_path: Path) -> Path:
    try:
        contract = json.loads(contract_path.read_text(encoding="utf-8-sig"))
    except (OSError, ValueError) as error:
        raise ValueError(f"物理設定JSONを読み取れません: {contract_path}") from error
    if not isinstance(contract, dict):
        raise ValueError("物理設定JSONはobject形式が必要です")
    phits_input = contract.get("phits_input")
    if not isinstance(phits_input, str) or not phits_input.strip():
        raise ValueError("物理設定JSONにphits_inputを指定してください")
    candidate = Path(phits_input).expanduser()
    if not candidate.is_absolute():
        candidate = contract_path.parent / candidate
    return _checked_input(candidate, "PHITS入力")


def _inspect_v2(path: Path) -> dict[str, object]:
    """Read only schema, counts, and provenance; never load deposit arrays."""
    import h5py

    with h5py.File(path, "r") as source:
        name = source.attrs.get("format", "")
        if isinstance(name, bytes):
            name = name.decode("utf-8")
        if name != "tes-dump2event" or int(source.attrs.get("format_version", -1)) != 2:
            raise ValueError("生成結果がHDF5 v2形式ではありません")
        for key in ("events", "particles", "tracks", "deposits", "provenance"):
            if key not in source:
                raise ValueError(f"HDF5 v2の/{key}がありません")
        events = source["events"]
        deposits = source["deposits"]
        count = len(events["event_id"])
        unresolved = 0
        # Event-level counters are small compared with the deposit ledger.
        for chunk_start in range(0, count, 10_000):
            unresolved += int(events["unresolved_count"][chunk_start:chunk_start + 10_000].sum())
        if unresolved:
            raise ValueError(f"未解決depositを含むファイルです: {unresolved}件")
        return {
            "events": count,
            "deposits": len(deposits["addition_id"]),
            "particles": len(source["particles"]["particle_id"]),
            "tracks": len(source["tracks"]["ncol"]),
            "phits_version": str(source["provenance"].attrs.get("phits_version", "unknown")),
            "algorithm": str(source["provenance"].attrs.get("algorithm_version", "unknown")),
            "resolved_sha256": str(source["provenance"].attrs.get("resolved_physics_sha256", "unavailable")),
            "unresolved": unresolved,
        }


def run_pb_egs5_ui(questionary, convert):
    """Prompt for real input paths, call the existing physics API, then inspect v2.

    The caller injects questionary and tes_cpp.dump2event so this module does
    not duplicate any physics and can be tested with small mock inputs.
    """
    print("PHITS 3.341 / Pb / EGS5 / 単位重み専用の実験的物理再構成です。")
    print("物理条件が異なる場合は変換が拒否されます。")
    print("HDF5 v2はPulse_Msへ直接渡せません。TES熱配分は未検証です。")

    dump_path = _ask_path(questionary, "無損失のDumpallファイルを指定してください")
    if dump_path is None:
        return None
    dump_path = _checked_input(dump_path, "Dumpall")

    config_path = _ask_path(questionary, "Pb物理設定JSONを指定してください")
    if config_path is None:
        return None
    config_path = _checked_input(config_path, "物理設定JSON")
    phits_input = _read_contract_input(config_path)
    print(f"設定に紐づくPHITS入力: {phits_input}")

    table_path = _ask_path(questionary, "privateなPb HATCH停止能CSVを指定してください")
    if table_path is None:
        return None
    table_path = _checked_input(table_path, "停止能CSV")

    output_path = _ask_path(
        questionary,
        "保存先HDF5 v2を指定してください（旧event.h5とは別名）",
        default=str(dump_path.parent / "event_pb_egs5_v2.h5"),
    )
    if output_path is None:
        return None
    if output_path.suffix.lower() not in (".h5", ".hdf5"):
        raise ValueError("出力ファイルは.h5または.hdf5を指定してください")
    if output_path.exists():
        raise FileExistsError(f"既存のファイルは上書きしません: {output_path}")
    if not output_path.parent.is_dir():
        raise FileNotFoundError(f"出力先ディレクトリが存在しません: {output_path.parent}")

    approved = questionary.confirm(
        "物理モードpb-egs5で変換しますか？（変換のみ。Pulse_Msは実行しません）",
        default=False,
    ).ask()
    if approved is not True:
        print("変換をキャンセルしました。")
        return None

    result = convert(
        dump_path,
        output_path,
        physics="pb-egs5",
        physics_config=config_path,
        stopping_table=table_path,
    )
    info = _inspect_v2(Path(result.output_path))
    print(
        f"HDF5 v2変換完了: {result.output_path} / "
        f"events={info['events']}, particles={info['particles']}, "
        f"tracks={info['tracks']}, deposits={info['deposits']}, "
        f"unresolved={info['unresolved']}"
    )
    print(f"provenance: PHITS={info['phits_version']}, algorithm={info['algorithm']}")
    print(f"resolved_physics_sha256={info['resolved_sha256']}")
    print("注意: スコア再構成の成功は全吸収やSEGMENT内の熱分配の証明ではありません。")
    return result
