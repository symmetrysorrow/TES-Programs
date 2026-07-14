import glob
import os
import re

def _list_child_dirs(path):
    return [entry.path for entry in os.scandir(path) if entry.is_dir()]


def _dirs_with_direct_dat(path):
    data_dirs = []
    for child in _list_child_dirs(path):
        if glob.glob(os.path.join(child, "*.dat")):
            data_dirs.append(child)
    return data_dirs


def _has_subdirs(path):
    return any(entry.is_dir() for entry in os.scandir(path))


def _has_pulse_layout(path):
    return any(
        entry.is_dir() and entry.name.endswith("_pulse")
        for entry in os.scandir(path)
    )


def _has_rt_layout(path):
    data_dirs = _dirs_with_direct_dat(path)
    if len(data_dirs) != 1:
        return False

    dats = glob.glob(os.path.join(data_dirs[0], "*.dat"))
    if not dats:
        return False
    sample = os.path.basename(dats[0])
    return re.search(r"^CH\d+_\d+mK_\d+uA\.dat$", sample) is not None


def _has_iv_layout(path):
    data_dirs = _dirs_with_direct_dat(path)
    if len(data_dirs) >= 2:
        return True
    if len(data_dirs) == 1 and _has_subdirs(data_dirs[0]):
        return True
    return False


def detect_dataset_kind(path):
    path = os.path.abspath(path)
    if not os.path.isdir(path):
        return None
    if _has_pulse_layout(path):
        return "pulse"
    if _has_iv_layout(path):
        return "iv"
    if _has_rt_layout(path):
        return "rt"
    return None


def dispatch(path, pulse_runner):
    kind = detect_dataset_kind(path)
    if kind == "pulse":
        pulse_runner(path)
        return kind
    if kind == "iv":
        from . import iv

        iv.main(path)
        return kind
    if kind == "rt":
        from . import rt as RT

        RT.main(path)
        return kind

    print("Could not detect the analysis target from the path.")
    print("Please choose a typical pulse, IV, or RT folder.")
    return None
