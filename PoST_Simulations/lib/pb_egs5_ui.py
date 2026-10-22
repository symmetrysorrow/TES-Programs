"""Interactive, opt-in PHITS Pb EGS5 -> HDF5 v2 conversion.

This UI intentionally does not call Pulse_Ms.  The v2 deposit schema has
SEGMENT records that require a separately validated TES allocation model.
"""

from __future__ import annotations

import json
import os
import hashlib
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
        incomplete = 0
        if "complete" not in events:
            raise ValueError("HDF5 v2の/events/completeがありません")
        # Read only small batches of event summaries; never load the deposits.
        for chunk_start in range(0, count, 10_000):
            sl = slice(chunk_start, chunk_start + 10_000)
            unresolved += int(events["unresolved_count"][sl].sum())
            complete = events["complete"][sl]
            incomplete += int(len(complete) - complete.sum())
        if unresolved or incomplete:
            raise ValueError(f"未解決または未完了のscoreがあります: unresolved={unresolved}, incomplete={incomplete}")
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


def _bundled_table_present() -> bool:
    """Report an installed private/bundled table only after its pinned hash matches."""
    name = "pb_egs5_phits3341.csv"
    expected = "20fd86eb0b57483bd19f20a8388d2c1ecb0431a02fa99f0933e2c3c4748514e9"
    repo_path = Path(__file__).resolve().parents[2] / "tes_cpp" / "dump2event" / "data" / name
    candidates = [repo_path]
    # scikit-build installs Python modules under <prefix>/python/tes_cpp and
    # native resources under <prefix>/tes_cpp/bin.
    candidates.extend(parent / "tes_cpp" / "bin" / name
                     for parent in Path(__file__).resolve().parents)
    try:
        from importlib.resources import files
    except (ImportError, ModuleNotFoundError):
        pass
    else:
        candidates.append(Path(str(files("tes_cpp").joinpath("bin", name))))
    # Environment configuration is an optional backward-compatible fallback.
    configured = os.environ.get("TES_CPP_PB_HATCH_TABLE")
    if configured:
        candidates.append(Path(configured).expanduser())
    for candidate in candidates:
        if candidate.is_file():
            digest = hashlib.sha256(candidate.read_bytes()).hexdigest()
            if digest != expected:
                raise ValueError(f"Pb HATCH CSV SHA-256 mismatch at {candidate}: {digest}")
            return True
    return False


def run_pb_egs5_ui(questionary, convert):
    """Prompt for real input paths, call the existing physics API, then inspect v2.

    The caller injects questionary and tes_cpp.dump2event so this module does
    not duplicate any physics and can be tested with small mock inputs.
    """
    print("PHITS 3.341 / Pb / EGS5物理再構成です。検証済み条件はPbのみです。")
    print("物理条件が異なる場合は変換が拒否されます。")
    print("HDF5 v2はPulse_Msへ直接渡せません。TES熱配分は未検証です。")

    dump_path = _ask_path(questionary, "無損失のDumpallファイルを指定してください")
    if dump_path is None:
        return None
    dump_path = _checked_input(dump_path, "Dumpall")

    print("Pb physics settings and the associated PHITS input are validated from the Dumpall run directory.")

    table_path = None
    if _bundled_table_present():
        print("検証済み同梱Pb HATCH停止能表を自動使用します（ドライバーでSHA-256照合）。")
    else:
        print("この環境には同梱Pb係数表がありません。元の非公開CSVを指定してください。")
        table_path = _ask_path(questionary, "元の検証済みPb HATCH停止能CSVを指定してください")
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
        **({"stopping_table": table_path} if table_path is not None else {}),
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
