import numpy as np
import argparse
import json
import pandas as pd
from pathlib import Path
from scipy.optimize import curve_fit
from scipy.interpolate import PchipInterpolator
from scipy.signal import bessel, filtfilt
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import re
import tqdm
import questionary
import h5py
from lib.pulse_hdf5 import iter_pulse_items as iter_hdf5_pulse_items


# Names used by the former directory-of-.dat output and the current JSON
# output.  The public target names are kept so existing calls still work.
PULSE_FILES = {
    "pulse_noise": "pulse_noise.json",
    "pulse_ms": "pulse_MS.json",
    "pulse_ms_noise": "pulse_MS_noise.json",
}

PULSE_HDF5_FILES = {
    "pulse_noise": "pulse_noise.h5",
    "pulse_ms": "pulse_MS.h5",
    "pulse_ms_noise": "pulse_MS_noise.h5",
}

# Every target accepted by the command line is also available in the
# interactive selector when --target is omitted.
TARGET_CHOICES = (
    "Pulse_ms",
    "Pulse_noise",
    "Pulse_ms_noise",
)
ALL_TARGETS_CHOICE = "All targets"
FULL_ENERGY_TARGETS = {"pulse_ms", "pulse_ms_noise"}
COMMON_MS_TARGETS = ("Pulse_ms", "Pulse_ms_noise")

# Pulse_ms position calibration is based on a short, symmetric average around
# the peak.  The old -10/+90 sample window mixed the fast CH0 decay and the
# slow CH1 rise, biasing CH1/CH0 toward the wrong absorber position.
PULSE_MS_PEAK_AVERAGE_HALF_WIDTH = 5


def pulse_file_name(target):
    """Return the current post_all JSON filename for *target*."""
    try:
        return PULSE_FILES[target.lower()]
    except KeyError as error:
        supported = ", ".join(sorted(PULSE_FILES))
        raise ValueError(f"Unknown pulse target {target!r}; supported targets: {supported}") from error


def pulse_data_path(data_path, position, target):
    """Return the available pulse file, preferring the current HDF5 format."""
    data_path = Path(data_path) / str(position)
    hdf5_path = data_path / PULSE_HDF5_FILES[target.lower()]
    if hdf5_path.is_file():
        return hdf5_path
    return data_path / pulse_file_name(target)


def output_csv_path(data_path, position, target, channel):
    """Feature CSV path in the target's dedicated numerical-results folder."""
    stem = Path(pulse_file_name(target)).stem
    return resolution_output_dir(data_path, target) / f"position_{position}" / f"{stem}_output_TES{channel}.csv"


def resolution_output_dir(data_path, target):
    """Directory for feature, resolution, and failure outputs of one target."""
    return Path(data_path) / "results" / target.lower()


def figure_output_dir(data_path, target):
    """Directory for plot images of one target."""
    return Path(data_path) / "figures" / target.lower()


def full_energy_event_ids(data_path, position, target, event_path=None):
    """Return the full-energy event IDs required for MS resolution analysis.

    ``FullEnergyList.dat`` is created alongside the external ``event.h5``
    data by ``PoST_Simulation.py``.  ``event_path`` is the external event
    root, matching ``PoST_Simulation.py``'s ``event_root``.
    """
    if target.lower() not in FULL_ENERGY_TARGETS:
        return None
    if event_path is None:
        raise ValueError(
            "EventPath is required for MS analysis. "
            "Specify the root containing position/<id>/FullEnergyList.dat."
        )
    event_root = Path(event_path).expanduser()
    path = event_root / str(position) / "FullEnergyList.dat"
    if not path.is_file():
        raise FileNotFoundError(
            f"Full-energy event list was not found for {target}: {path}. "
            "Check EventPath and run Dump2Event before extracting MS features."
        )
    event_ids = {
        line.strip()
        for line in path.read_text(encoding="utf-8").splitlines()
        if line.strip()
    }
    if not event_ids:
        raise ValueError(
            f"Full-energy event list is empty: {path}. "
            "Check that EventPath matches the simulation energy and rerun "
            "Dump2Event with the correct input.json."
        )
    return event_ids


def bessel_filter(pulse, rate, cutoff):
    """Apply the noise-target low-pass filter without depending on getpara."""
    normalized_cutoff = float(cutoff) / float(rate) * 2.0
    coefficients_b, coefficients_a = bessel(2, normalized_cutoff, "low")
    return filtfilt(coefficients_b, coefficients_a, pulse)


class JsonStream:
    """Incremental JSON reader which does not retain skipped values in memory."""

    def __init__(self, path, chunk_size=1024 * 1024):
        self.file = open(path, "r", encoding="utf-8")
        self.chunk_size = chunk_size
        self.buffer = ""
        self.position = 0
        self.eof = False
        self.decoder = json.JSONDecoder()

    def close(self):
        self.file.close()

    def _read_more(self):
        # Discard consumed input before adding a chunk.  This is essential for
        # multi-gigabyte ``time`` arrays and pulse documents.
        if self.position:
            self.buffer = self.buffer[self.position:]
            self.position = 0
        chunk = self.file.read(self.chunk_size)
        if chunk:
            self.buffer += chunk
        else:
            self.eof = True

    def _ensure_character(self):
        while self.position >= len(self.buffer) and not self.eof:
            self._read_more()
        if self.position >= len(self.buffer):
            raise ValueError("Unexpected end of JSON input")

    def peek(self):
        self._ensure_character()
        return self.buffer[self.position]

    def get(self):
        character = self.peek()
        self.position += 1
        return character

    def whitespace(self):
        while not self.eof or self.position < len(self.buffer):
            try:
                if not self.peek().isspace():
                    return
                self.position += 1
            except ValueError:
                return

    def expect(self, expected):
        self.whitespace()
        actual = self.get()
        if actual != expected:
            raise ValueError(f"Expected {expected!r}, found {actual!r}")

    def decode_value(self):
        """Decode one value, retaining only that value's JSON text."""
        self.whitespace()
        while True:
            try:
                value, end = self.decoder.raw_decode(self.buffer, self.position)
                self.position = end
                return value
            except json.JSONDecodeError as error:
                if self.eof:
                    raise ValueError(f"Invalid or truncated JSON near offset {self.position}") from error
                self._read_more()

    def skip_value(self):
        """Consume one JSON value without constructing a Python object."""
        self.whitespace()
        first = self.peek()
        if first == '"':
            self.decode_value()
            return
        if first not in "[{":
            while True:
                character = self.peek()
                if character.isspace() or character in ",}]":
                    return
                self.position += 1

        stack = [self.get()]
        in_string = False
        escaped = False
        while stack:
            character = self.get()
            if in_string:
                if escaped:
                    escaped = False
                elif character == "\\":
                    escaped = True
                elif character == '"':
                    in_string = False
            elif character == '"':
                in_string = True
            elif character in "[{":
                stack.append(character)
            elif character in "]}":
                opener = stack.pop()
                if (opener, character) not in (("[", "]"), ("{", "}")):
                    raise ValueError("Mismatched JSON brackets")


def iter_pulse_items(path):
    """Yield ``(event_id, pulse)`` from post_all JSON without loading it all.

    ``input`` and the potentially enormous shared ``time`` array are skipped
    at byte-stream level.  Only one event's CH0/CH1 arrays exists in memory.
    """
    stream = JsonStream(path)
    found_pulses = False
    event_count = 0
    try:
        stream.expect("{")
        while True:
            stream.whitespace()
            if stream.peek() == "}":
                stream.get()
                break
            key = stream.decode_value()
            if not isinstance(key, str):
                raise ValueError(f"{path}: JSON object key is not a string")
            stream.expect(":")
            if key == "pulses":
                found_pulses = True
                stream.expect("{")
                while True:
                    stream.whitespace()
                    if stream.peek() == "}":
                        stream.get()
                        break
                    event_id = stream.decode_value()
                    if not isinstance(event_id, str):
                        raise ValueError(f"{path}: pulse ID is not a string")
                    stream.expect(":")
                    yield event_id, stream.decode_value()
                    event_count += 1
                    stream.whitespace()
                    separator = stream.get()
                    if separator == "}":
                        break
                    if separator != ",":
                        raise ValueError(f"{path}: expected ',' or '}}' in pulses")
            else:
                stream.skip_value()

            stream.whitespace()
            separator = stream.get()
            if separator == "}":
                break
            if separator != ",":
                raise ValueError(f"{path}: expected ',' or '}}' in root object")
    finally:
        stream.close()

    if not found_pulses:
        raise ValueError(f"{path} must contain a 'pulses' object")
    if event_count == 0:
        raise ValueError(f"{path} contains no pulses")


def save_readpulse_debug(Data_path, target, pulse_raw, pulse_filt, para, path, reason, peak_index=None, rise_10=None, rise_90=None):
    debug_dir = figure_output_dir(Data_path, target) / "debug_readpulse"
    debug_dir.mkdir(parents=True, exist_ok=True)

    t = np.arange(len(pulse_raw)) / para["rate"]
    st_l = max(0, para["SettlingTime"] - 10)
    st_r = min(len(pulse_filt), para["SettlingTime"] + 90)

    plt.figure(figsize=(10, 5))
    plt.plot(t, pulse_raw, color="gray", alpha=0.35, label="raw")
    plt.plot(t, pulse_filt, color="navy", linewidth=1.2, label="filtered")

    if peak_index is not None and 0 <= peak_index < len(t):
        plt.axvline(t[peak_index], color="red", linestyle="--", label="peak")
    if rise_10 is not None and 0 <= rise_10 < len(t):
        plt.axvline(t[rise_10], color="green", linestyle="--", label="10%")
    if rise_90 is not None and 0 <= rise_90 < len(t):
        plt.axvline(t[rise_90], color="orange", linestyle="--", label="90%")
    if len(pulse_filt) > 0:
        peak = np.nanmax(pulse_filt)
        plt.axhline(peak * 0.1, color="green", alpha=0.25)
        plt.axhline(peak * 0.9, color="orange", alpha=0.25)
    plt.axvspan(st_l / para["rate"], st_r / para["rate"], color="gray", alpha=0.18, label="Settling window")
    plt.title(f"{path}\nreason: {reason}")
    plt.xlabel("Time [s]")
    plt.ylabel("Amplitude")
    plt.grid(True)
    plt.legend(fontsize=8)
    plt.tight_layout()

    safe_path = re.sub(r'[\\/:*?"<>|]', "_", str(path))
    plt.savefig(debug_dir / f"{safe_path}.png", dpi=200)
    plt.close()


def robust_scale(values):
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if len(values) == 0:
        return np.nan, np.nan

    median = np.median(values)
    mad = np.median(np.abs(values - median))
    if np.isfinite(mad) and mad > 0:
        scale = 1.4826 * mad
    else:
        q25, q75 = np.percentile(values, [25, 75])
        iqr = q75 - q25
        scale = iqr / 1.349 if np.isfinite(iqr) and iqr > 0 else np.std(values)

    return median, scale


def remove_outlier_events(CH0_df, CH1_df, z_thresh=4.5):
    """
    Remove event pairs whose CH0 pulse shape is clearly off.

    CH0 is used as the primary quality gate here because the abnormal
    17-position behavior is concentrated on CH0. We keep the two channel
    CSV rows aligned by dropping the same event indices from both.
    """
    mask = (
        np.isfinite(CH0_df["height"].to_numpy(dtype=float))
        & np.isfinite(CH0_df["peak_index"].to_numpy(dtype=float))
        & np.isfinite(CH1_df["height"].to_numpy(dtype=float))
    )

    energy_column = "energy_height" if "energy_height" in CH0_df and "energy_height" in CH1_df else "height"
    if energy_column != "height":
        mask &= (
            np.isfinite(CH0_df[energy_column].to_numpy(dtype=float))
            & np.isfinite(CH1_df[energy_column].to_numpy(dtype=float))
        )
    sum_height = CH0_df[energy_column].to_numpy(dtype=float) + CH1_df[energy_column].to_numpy(dtype=float)
    log_sum_height = np.log10(np.clip(sum_height, np.finfo(float).tiny, None))

    median, scale = robust_scale(log_sum_height[mask])
    if np.isfinite(scale) and scale > 0:
        z = np.abs(log_sum_height - median) / scale
        mask &= np.isfinite(z) & (z <= z_thresh)

    median, scale = robust_scale(CH0_df.loc[mask, "peak_index"].to_numpy(dtype=float))
    if np.isfinite(scale) and scale > 0:
        z = np.abs(CH0_df["peak_index"].to_numpy(dtype=float) - median) / scale
        mask &= np.isfinite(z) & (z <= z_thresh)

    return CH0_df.loc[mask].reset_index(drop=True), CH1_df.loc[mask].reset_index(drop=True), mask


def ReadPulse(
    Data_path,
    pulse,
    target="Pulse_ms",
    path="",
    debug_plot=False,
    failures=None,
    processing_filter=None,
    peak_average_half_width=None,
    settling_half_width=None,
    fixed_peak_index=None,
):
    with open(f"{Data_path}/input.json") as f:
        para = json.load(f)

    # Callers often retain the original waveform (e.g. CH0/CH1 in one JSON
    # event), so never change it in place while normalising its polarity.
    pulse = np.asarray(pulse, dtype=float).copy()
    if np.mean(pulse) <= 0:
        pulse *= -1

    if processing_filter is None:
        use_bessel = target.lower() != "pulse_ms"
    elif processing_filter.lower() in {"bessel", "bessel_filter"}:
        use_bessel = True
    elif processing_filter.lower() in {"none", "raw"}:
        use_bessel = False
    else:
        raise ValueError(f"Unknown processing filter: {processing_filter!r}")

    pulse_filt = bessel_filter(pulse, para["rate"], para["cutoff"]) if use_bessel else pulse

    def failed(reason):
        if failures is not None:
            failures.append({"event": str(path).rsplit(":", 1)[-1], "reason": reason})
        if debug_plot:
            save_readpulse_debug(Data_path, target, pulse, pulse_filt, para, path, reason)
        return [np.nan, np.nan, np.nan, np.nan, np.nan]

    if not np.all(np.isfinite(pulse_filt)):
        return failed("non-finite values")

    peak = np.max(pulse_filt)
    peak_index = int(np.argmax(pulse_filt))

    if peak_index < 10:
        return failed("peak too early")

    if peak_average_half_width is not None:
        half_width = int(peak_average_half_width)
        if half_width < 0:
            raise ValueError("peak_average_half_width must be non-negative")
        l = max(0, peak_index - half_width)
        r = min(len(pulse_filt), peak_index + half_width + 1)
    elif target.lower() == "pulse_ms":
        half_width = PULSE_MS_PEAK_AVERAGE_HALF_WIDTH
        l = max(0, peak_index - half_width)
        r = min(len(pulse_filt), peak_index + half_width + 1)
    else:
        # Keep the existing feature definition for the separately filtered
        # noise target until it has its own matching calibration curve.
        l = max(0, peak_index - 10)
        r = min(len(pulse_filt), peak_index + 90)
    peak_av = np.mean(pulse_filt[l:r])

    if fixed_peak_index is None:
        energy_height = peak_av
    else:
        if isinstance(fixed_peak_index, (bool, np.bool_)):
            return failed("fixed_peak_index must be an integer")
        try:
            fixed_peak_index_float = float(fixed_peak_index)
        except (TypeError, ValueError):
            return failed("fixed_peak_index must be an integer")
        if not np.isfinite(fixed_peak_index_float) or fixed_peak_index_float != int(fixed_peak_index_float):
            return failed("fixed_peak_index must be an integer")
        fixed_peak_index = int(fixed_peak_index_float)
        if fixed_peak_index < 0 or fixed_peak_index >= len(pulse_filt):
            return failed("fixed_peak_index outside waveform")
        half_width = PULSE_MS_PEAK_AVERAGE_HALF_WIDTH
        fixed_l = max(0, fixed_peak_index - half_width)
        fixed_r = min(len(pulse_filt), fixed_peak_index + half_width + 1)
        if fixed_l >= fixed_r:
            return failed("invalid fixed energy window")
        energy_height = np.mean(pulse_filt[fixed_l:fixed_r])

    rise_90 = None
    for i in reversed(range(0, peak_index)):
        if pulse_filt[i] <= peak * 0.9:
            rise_90 = i
            break
    if rise_90 is None:
        return failed("rise_90 not found")

    rise_10 = None
    for j in reversed(range(0, rise_90)):
        if pulse_filt[j] <= peak * 0.1:
            rise_10 = j
            break
    if rise_10 is None:
        return failed("rise_10 not found")

    rise = (rise_90 - rise_10) / para["rate"]

    if settling_half_width is None:
        st_l = para["SettlingTime"] - 10
        st_r = para["SettlingTime"] + 90
    else:
        settling_half_width = int(settling_half_width)
        if settling_half_width < 0:
            raise ValueError("settling_half_width must be non-negative")
        st_l = max(0, peak_index - settling_half_width)
        st_r = min(len(pulse_filt), peak_index + settling_half_width + 1)
    if st_l < 0 or st_r > len(pulse_filt) or st_l >= st_r:
        return failed("invalid settling window")

    ST_window = pulse_filt[st_l:st_r]
    if len(ST_window) == 0 or not np.all(np.isfinite(ST_window)):
        return failed("invalid settling samples")

    ST_height = np.mean(ST_window)

    if debug_plot:
        t = np.arange(len(pulse_filt)) / para["rate"]
        plt.figure(figsize=(8, 4))
        plt.plot(t, pulse_filt, label="pulse")
        plt.axvline(t[peak_index], color="r", linestyle="--", label="peak")
        plt.axvline(t[rise_10], color="g", linestyle="--", label="10%")
        plt.axvline(t[rise_90], color="orange", linestyle="--", label="90%")
        plt.axhline(peak * 0.1, color="g", alpha=0.3)
        plt.axhline(peak * 0.9, color="orange", alpha=0.3)
        plt.axvspan(st_l / para["rate"], st_r / para["rate"], color="gray", alpha=0.2, label="Settling window")
        plt.xlabel("Time [s]")
        plt.ylabel("Amplitude")
        plt.legend()
        plt.tight_layout()
        plt.show()

    return [peak_av, peak_index, rise, ST_height, energy_height]


def reference_peak_indices(data_path, position, target="Pulse_ms", processing_filter=None):
    """Return filtered noiseless-reference peak indices for CH0 and CH1."""
    path = Path(data_path) / "pulses.h5"
    if not path.is_file():
        raise FileNotFoundError(f"Reference pulse file was not found: {path}")
    try:
        reference = next(iter_hdf5_pulse_items(path, event_ids={str(position)}))[1]
    except StopIteration as error:
        raise ValueError(f"Reference pulses.h5 does not contain position {position}") from error
    with open(Path(data_path) / "input.json", encoding="utf-8") as file:
        para = json.load(file)
    if processing_filter is None:
        processing_filter = "none" if target.lower() == "pulse_ms" else "bessel"
    indices = []
    for channel in ("ch0", "ch1"):
        pulse = np.asarray(reference[channel], dtype=float).copy()
        if np.mean(pulse) <= 0:
            pulse *= -1
        if processing_filter.lower() in {"bessel", "bessel_filter"}:
            pulse = bessel_filter(pulse, para["rate"], para["cutoff"])
        elif processing_filter.lower() not in {"none", "raw"}:
            raise ValueError(f"Unknown processing filter: {processing_filter!r}")
        if pulse.size == 0 or not np.all(np.isfinite(pulse)):
            raise ValueError(f"Invalid reference waveform for position {position}, {channel}")
        indices.append(int(np.argmax(pulse)))
    return tuple(indices)


# Descriptive alias for callers that want to make the reference role explicit.
load_reference_peak_indices = reference_peak_indices


def MakeOutput(Data_path, target, event_path=None):
    """Extract pulse features from the current post_all HDF5 or JSON files.

    HDF5 is preferred when available, with JSON retained as a fallback.
    MS targets retain only FullEnergyList.dat events before feature extraction.
    Feature CSVs are written to the target's numerical-results folder.
    """
    data_path = Path(Data_path)
    with open(data_path / "input.json", encoding="utf-8") as f:
        para = json.load(f)

    if target.lower() in FULL_ENERGY_TARGETS and event_path is None:
        event_path = ask_event_path(data_path)
        if event_path is None:
            raise SystemExit(0)

    if target.lower() == "pulse_ms":
        print("Pulse_ms: no Bessel filter")
    else:
        print("BesselFilter")

    result_dir = resolution_output_dir(data_path, target)
    result_dir.mkdir(parents=True, exist_ok=True)
    feature_summary = []

    for i, posi in enumerate(para["position"]):
        print(f"{i+1}/{len(para['position'])}")
        source_path = pulse_data_path(data_path, posi, target)
        if not source_path.is_file():
            raise FileNotFoundError(f"Pulse data was not found: {source_path}")
        full_energy_ids = full_energy_event_ids(data_path, posi, target, event_path)
        fixed_indices = reference_peak_indices(data_path, posi, target)
        results = {0: [], 1: []}
        pulse_ids = []
        failures = []
        source_event_count = None
        total = None
        if source_path.suffix.lower() == ".h5":
            with h5py.File(source_path, "r") as file:
                source_event_count = len(file["event_id"])
                total = source_event_count
        pulse_iterator = (
            iter_hdf5_pulse_items(source_path, event_ids=full_energy_ids)
            if source_path.suffix.lower() == ".h5"
            else iter_pulse_items(source_path)
        )
        for pulse_id, pulse_data in tqdm.tqdm(
            pulse_iterator, total=total, desc=f"Position {posi}", unit="event"
        ):
            if full_energy_ids is not None and str(pulse_id) not in full_energy_ids:
                continue
            if not isinstance(pulse_data, dict):
                raise ValueError(f"Pulse {pulse_id!r} in {source_path} is not an object")
            pulse_ids.append(pulse_id)
            for ch in [0, 1]:
                try:
                    pulse = pulse_data[f"ch{ch}"]
                except KeyError as error:
                    raise ValueError(
                        f"Pulse {pulse_id!r} in {source_path} has no ch{ch} waveform"
                    ) from error
                results[ch].append(
                    ReadPulse(
                        Data_path, pulse, path=f"{source_path}:{pulse_id}",
                        target=target, failures=failures,
                        fixed_peak_index=fixed_indices[ch],
                    )
                )

        columns = ["height", "peak_index", "rise", "ST_Height", "energy_height"]
        feature_frames = {}
        for ch in [0, 1]:
            df = pd.DataFrame(results[ch], columns=columns, index=pulse_ids)
            df.index.name = "id"
            feature_frames[ch] = df
            output_path = output_csv_path(data_path, posi, target, ch)
            output_path.parent.mkdir(parents=True, exist_ok=True)
            df.to_csv(output_path)
        failure_path = resolution_output_dir(data_path, target) / f"reso_failures_{target}_position_{posi}.json"
        failure_path.parent.mkdir(parents=True, exist_ok=True)
        with open(failure_path, "w", encoding="utf-8") as file:
            json.dump(failures, file, ensure_ascii=False, indent=2)
        if failures:
            failed_events = {item["event"] for item in failures}
            print(f"Position {posi}: Failed {len(failed_events)} events ({len(failures)} channel analyses)")
        else:
            failed_events = set()
        valid_feature_events = (
            np.isfinite(feature_frames[0]["height"])
            & np.isfinite(feature_frames[0]["peak_index"])
            & np.isfinite(feature_frames[1]["height"])
            & np.isfinite(feature_frames[0]["energy_height"])
            & np.isfinite(feature_frames[1]["energy_height"])
        ).sum()
        feature_summary.append({
            "position": posi,
            "source_event_count": source_event_count if source_event_count is not None else len(pulse_ids),
            "selected_event_count": len(pulse_ids),
            "full_energy_event_count": len(full_energy_ids) if full_energy_ids is not None else np.nan,
            "failed_event_count": len(failed_events),
            "valid_feature_event_count": int(valid_feature_events),
        })

    pd.DataFrame(feature_summary).to_csv(
        result_dir / f"feature_summary_{target}.csv", index=False
    )


def gaussian(x, amp, mean, stddev):
    return amp * np.exp(-((x - mean) ** 2) / (2 * stddev ** 2))


def optimal_bin_count(data):
    q1, q3 = np.percentile(data, [25, 75])
    iqr = q3 - q1
    if iqr <= 0:
        return 40
    bin_width = 2 * iqr / (len(data) ** (1 / 3))
    if not np.isfinite(bin_width) or bin_width <= 0:
        return 40
    bin_count = int(np.ceil((np.max(data) - np.min(data)) / bin_width))
    return int(np.clip(bin_count, 20, 40))


def MakeHistgram(
    data,
    posi,
    HistColor=None,
    bin_num=None,
    stabilize=True,
    fit_method="histogram",
):
    """Plot a histogram and return FWHM and relative resolution.

    ``fit_method='unbinned_gaussian'`` uses the Gaussian maximum-likelihood
    estimates of the finite samples.  The histogram remains a display only,
    so the result does not depend on bin edges or bin count.
    """
    data = np.asarray(data, dtype=float)
    data = data[np.isfinite(data)]
    if len(data) < 3 or np.ptp(data) == 0:
        return np.nan, np.nan

    fit_method = fit_method.lower()
    bin_num = optimal_bin_count(data) if bin_num is None else int(bin_num)
    hist, bin_edges = np.histogram(data, bins=bin_num, density=False)
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2

    if HistColor is not None:
        plt.hist(data, bins=bin_num, density=False, label=f"abs-{posi}", color=HistColor)
    else:
        plt.hist(data, bins=bin_num, density=False, label=f"abs-{posi}")

    if fit_method in {"unbinned", "unbinned_gaussian"}:
        mean_fit = np.mean(data)
        stddev_fit = np.std(data, ddof=0)
        if not np.isfinite(mean_fit) or not np.isfinite(stddev_fit) or stddev_fit <= 0:
            return np.nan, np.nan
        fwhm = 2 * stddev_fit * np.sqrt(2 * np.log(2))
        bin_width = np.mean(np.diff(bin_edges))
        amplitude = len(data) * bin_width / (stddev_fit * np.sqrt(2 * np.pi))
        x_fit = np.linspace(bin_edges[0], bin_edges[-1], 1000)
        plt.plot(x_fit, gaussian(x_fit, amplitude, mean_fit, stddev_fit), color="red", alpha=0.5)
        return fwhm, fwhm / mean_fit if mean_fit > 0 else np.nan

    if fit_method not in {"histogram", "binned"}:
        raise ValueError(f"Unknown histogram fit method: {fit_method!r}")

    initial_guess = [np.max(hist), np.mean(data), np.std(data)]

    try:
        popt, _ = curve_fit(gaussian, bin_centers, hist, p0=initial_guess, maxfev=100000)
    except (RuntimeError, ValueError, FloatingPointError):
        # A histogram fit can be ill-conditioned for a sparse distribution.
        # Report a direct width instead of aborting every position's analysis.
        mean_fit = np.mean(data)
        fwhm = 2 * np.std(data, ddof=1) * np.sqrt(2 * np.log(2))
        return fwhm, fwhm / mean_fit if mean_fit > 0 else np.nan
    _, mean_fit, stddev_fit = popt
    fwhm = 2 * stddev_fit * np.sqrt(2 * np.log(2))

    reference_width = 2 * np.std(data) * np.sqrt(2 * np.log(2))
    sample_mean = np.mean(data)
    robust_sigma = (np.percentile(data, 84) - np.percentile(data, 16)) / 2
    robust_fwhm = 2 * robust_sigma * np.sqrt(2 * np.log(2))
    robust_reso = robust_fwhm / sample_mean if sample_mean > 0 else np.nan
    fit_reso = fwhm / mean_fit if mean_fit > 0 else np.nan

    if stabilize and np.isfinite(reference_width) and reference_width > 0 and fwhm < 0.6 * reference_width:
        retry_bins = min(max(12, bin_num // 2), 14)
        if retry_bins != bin_num:
            hist, bin_edges = np.histogram(data, bins=retry_bins, density=False)
            bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
            initial_guess = [np.max(hist), np.mean(data), np.std(data)]
            try:
                popt, _ = curve_fit(gaussian, bin_centers, hist, p0=initial_guess, maxfev=100000)
                _, mean_fit, stddev_fit = popt
                fwhm = 2 * stddev_fit * np.sqrt(2 * np.log(2))
                fit_reso = fwhm / mean_fit if mean_fit > 0 else np.nan
            except (RuntimeError, ValueError, FloatingPointError):
                fwhm = robust_fwhm
                mean_fit = sample_mean

    # If the Gaussian fit becomes much broader than a robust percentile-based
    # estimate, fall back to the robust width so one skewed histogram binning
    # does not dominate the reported resolution.
    if stabilize and np.isfinite(robust_reso) and np.isfinite(fit_reso) and fit_reso > 1.5 * robust_reso:
        fwhm = robust_fwhm
        mean_fit = sample_mean

    x_fit = np.linspace(bin_edges[0], bin_edges[-1], 1000)
    plt.plot(x_fit, gaussian(x_fit, *popt), color="red", alpha=0.5)

    return fwhm, fwhm / mean_fit


def generate_symmetric_colors(n):
    half_n = n // 2
    hues = np.linspace(0.1, 0.9, half_n)
    hues = np.concatenate([hues, hues[::-1]])
    if n % 2 != 0:
        hues = np.concatenate([hues[:half_n], [0.0], hues[half_n:]])
    colors = [mcolors.hsv_to_rgb((h, 1.0, 1.0)) for h in hues]
    return np.array(colors)


def load_ratio_calibration(data_path):
    """Return a bounded interpolator from current ``ratios.csv`` data.

    New files use the same short peak-average definition as Pulse_ms.  The
    old maximum-ratio column remains accepted for older result directories.
    """
    ratio_path = Path(data_path) / "ratios.csv"
    if not ratio_path.is_file():
        raise FileNotFoundError(f"Ratio calibration was not found: {ratio_path}")

    table = pd.read_csv(ratio_path)
    ratio_column = (
        "ch1_ch0_peak_average_ratio"
        if "ch1_ch0_peak_average_ratio" in table.columns
        else "ch1_ch0_max_ratio"
    )
    required = {"x_mm", ratio_column}
    missing = required - set(table.columns)
    if missing:
        raise ValueError(f"{ratio_path} is missing columns: {', '.join(sorted(missing))}")
    table = table.dropna(subset=["x_mm", ratio_column]).sort_values(ratio_column)
    table = table.drop_duplicates(ratio_column)
    if len(table) < 2:
        raise ValueError(f"{ratio_path} needs at least two distinct ratio samples")

    ratios = table[ratio_column].to_numpy(float)
    positions = table["x_mm"].to_numpy(float)
    return PchipInterpolator(ratios, positions, extrapolate=False)


def load_feature_pair(data_path, position, target, event_path=None):
    """Load matching channel features and retain only shared event IDs."""
    ch0_path = output_csv_path(data_path, position, target, 0)
    ch1_path = output_csv_path(data_path, position, target, 1)
    if not ch0_path.is_file() or not ch1_path.is_file():
        raise FileNotFoundError(
            f"Feature CSVs were not found for position {position}. Run MakeOutput first."
        )
    ch0 = pd.read_csv(ch0_path, index_col="id")
    ch1 = pd.read_csv(ch1_path, index_col="id")
    common_ids = ch0.index.intersection(ch1.index, sort=False)
    if common_ids.empty:
        raise ValueError(f"No matching event IDs in {ch0_path} and {ch1_path}")
    if target.lower() in FULL_ENERGY_TARGETS:
        allowed_ids = full_energy_event_ids(data_path, position, target, event_path)
        non_full_ids = common_ids[~common_ids.astype(str).isin(allowed_ids)]
        if len(non_full_ids):
            raise ValueError(
                f"Feature CSVs for {target}, position {position} include "
                f"{len(non_full_ids)} non-full-energy events. Run MakeOutput again."
            )
    filtered_ch0, filtered_ch1, mask = remove_outlier_events(
        ch0.loc[common_ids], ch1.loc[common_ids]
    )
    return filtered_ch0, filtered_ch1, common_ids[~mask].astype(str).tolist()


def _read_common_target_features(data_path, position, target, event_path, event_ids):
    """Read one target using the processing settings of the fair MS comparison."""
    source_path = pulse_data_path(data_path, position, target)
    if not source_path.is_file():
        raise FileNotFoundError(f"Pulse data was not found: {source_path}")

    if source_path.suffix.lower() == ".h5":
        pulse_iterator = iter_hdf5_pulse_items(source_path, event_ids=event_ids)
        with h5py.File(source_path, "r") as file:
            total = len(event_ids) if event_ids is not None else len(file["event_id"])
    else:
        pulse_iterator = iter_pulse_items(source_path)
        total = len(event_ids) if event_ids is not None else None

    features = {}
    fixed_indices = reference_peak_indices(data_path, position, target, processing_filter="bessel")
    progress = tqdm.tqdm(
        pulse_iterator,
        total=total,
        desc=f"Position {position}: {target}",
        position=1,
        leave=False,
        dynamic_ncols=True,
    )
    for pulse_id, pulse_data in progress:
        pulse_id = str(pulse_id)
        if event_ids is not None and pulse_id not in event_ids:
            continue
        values = []
        for channel in (0, 1):
            values.append(
                ReadPulse(
                    data_path,
                    pulse_data[f"ch{channel}"],
                    target=target,
                    path=f"{source_path}:{pulse_id}",
                    # Use the common Bessel filter. Position ratios retain
                    # adaptive 11-sample peaks; energy uses fixed 11-sample
                    # windows at the noiseless-reference peak.
                    processing_filter="bessel",
                    peak_average_half_width=PULSE_MS_PEAK_AVERAGE_HALF_WIDTH,
                    fixed_peak_index=fixed_indices[channel],
                )
            )
        features[pulse_id] = values
    return features


def CommonMSComparison(Data_path, position, event_path, show=False, bin_num=None):
    """Compare MS and MS+noise using common event/filter conventions.

    Position ratios retain adaptive 11-sample peaks. Energy Max/Sum/Min use
    fixed 11-sample noiseless-reference windows; ST uses SettlingTime.

    ``Pulse_noise`` is reported as an additional, independent baseline.  It
    is a separate single-site event population, so it has no event-by-event
    counterpart to MS and is not part of the common-event comparison.
    Results are written separately from the ordinary target-specific ``Resos``
    output so the latter remains backward compatible.
    """
    data_path = Path(Data_path)
    position = int(position)
    with open(data_path / "input.json", encoding="utf-8") as file:
        para = json.load(file)
    if position not in para["position"]:
        raise ValueError(f"Position {position} is not in input.json: {para['position']}")
    if event_path is None:
        raise ValueError("EventPath is required for the common MS comparison")

    event_ids = full_energy_event_ids(data_path, position, "Pulse_ms", event_path)
    target_features = {
        target: _read_common_target_features(
            data_path, position, target, event_path, event_ids
        )
        for target in COMMON_MS_TARGETS
    }
    common_ids = [
        event_id
        for event_id in target_features[COMMON_MS_TARGETS[0]]
        if event_id in target_features[COMMON_MS_TARGETS[1]]
    ]
    if not common_ids:
        raise ValueError(f"No common FullEnergy events found at position {position}")

    valid_ids = []
    estimators = ("Max", "Sum", "Min", "ST")
    values = {
        target: {estimator: [] for estimator in estimators}
        for target in COMMON_MS_TARGETS
    }
    for event_id in common_ids:
        row = {
            target: target_features[target][event_id]
            for target in COMMON_MS_TARGETS
        }
        if not all(np.all(np.isfinite(row[target][channel][0:5]))
                   for target in COMMON_MS_TARGETS for channel in (0, 1)):
            continue
        valid_ids.append(event_id)
        for target in COMMON_MS_TARGETS:
            ch0, ch1 = row[target]
            values[target]["Max"].append(max(ch0[4], ch1[4]))
            values[target]["Sum"].append(ch0[4] + ch1[4])
            values[target]["Min"].append(min(ch0[4], ch1[4]))
            values[target]["ST"].append(ch0[3] + ch1[3])

    if len(valid_ids) < 3:
        raise ValueError(
            f"Fewer than three valid common events at position {position}: "
            f"{len(valid_ids)}"
        )

    position_from_ratio = load_ratio_calibration(data_path)
    position_values = {target: [] for target in COMMON_MS_TARGETS}
    for event_id in valid_ids:
        for target in COMMON_MS_TARGETS:
            ch0, ch1 = target_features[target][event_id]
            ratio = ch1[0] / ch0[0] if ch0[0] != 0 else np.nan
            position_values[target].append(position_from_ratio(ratio))

    noise_values = None
    noise_source = pulse_data_path(data_path, position, "Pulse_noise")
    if noise_source.is_file():
        noise_features = _read_common_target_features(
            data_path, position, "Pulse_noise", event_path, None
        )
        noise_values = {estimator: [] for estimator in estimators}
        for channel_features in noise_features.values():
            if not all(np.all(np.isfinite(channel_features[channel][0:5])) for channel in (0, 1)):
                continue
            ch0, ch1 = channel_features
            noise_values["Max"].append(max(ch0[4], ch1[4]))
            noise_values["Sum"].append(ch0[4] + ch1[4])
            noise_values["Min"].append(min(ch0[4], ch1[4]))
            noise_values["ST"].append(ch0[3] + ch1[3])

    output_dir = data_path / "results" / "common_comparison"
    figure_dir = data_path / "figures" / "common_comparison"
    output_dir.mkdir(parents=True, exist_ok=True)
    figure_dir.mkdir(parents=True, exist_ok=True)

    rows = []
    for estimator in estimators:
        for target in COMMON_MS_TARGETS:
            target_label = target.lower()
            plt.figure()
            fwhm_current, resolution_fraction = MakeHistgram(
                values[target][estimator],
                position,
                bin_num=bin_num,
                stabilize=False,
                fit_method="unbinned_gaussian",
            )
            plt.title(f"{estimator}: position {position}")
            save_histogram(
                figure_dir / f"position_{position}_{estimator.lower()}_{target_label}.png",
                "Current [A]",
                f"{estimator}: position {position}",
                show,
            )
            array = np.asarray(values[target][estimator], dtype=float)
            rows.append({
                "position": position,
                "metric": "energy",
                "estimator": estimator,
                "target": target,
                "event_set": "common_full_energy",
                "fit_method": "unbinned_gaussian",
                "event_count": len(valid_ids),
                "fwhm_current": fwhm_current,
                "resolution_fraction": resolution_fraction,
                "resolution_keV": resolution_fraction * float(para["E"]),
                "mean_current": np.mean(array),
                "std_current": np.std(array, ddof=1),
            })

    if noise_values is not None:
        for estimator in estimators:
            if len(noise_values[estimator]) < 3:
                continue
            plt.figure()
            fwhm_current, resolution_fraction = MakeHistgram(
                noise_values[estimator],
                position,
                bin_num=bin_num,
                stabilize=False,
                fit_method="unbinned_gaussian",
            )
            plt.title(f"{estimator}: position {position}")
            save_histogram(
                figure_dir / f"position_{position}_{estimator.lower()}_pulse_noise.png",
                "Current [A]",
                f"{estimator}: position {position}",
                show,
            )
            array = np.asarray(noise_values[estimator], dtype=float)
            rows.append({
                "position": position,
                "metric": "energy",
                "estimator": estimator,
                "target": "Pulse_noise",
                "event_set": "independent_noise",
                "fit_method": "unbinned_gaussian",
                "event_count": len(array),
                "fwhm_current": fwhm_current,
                "resolution_fraction": resolution_fraction,
                "resolution_keV": resolution_fraction * float(para["E"]),
                "mean_current": np.mean(array),
                "std_current": np.std(array, ddof=1),
            })

    for target in COMMON_MS_TARGETS:
        target_label = target.lower()
        plt.figure()
        position_fwhm, _ = MakeHistgram(
            position_values[target],
            position,
            bin_num=bin_num,
            stabilize=False,
            fit_method="unbinned_gaussian",
        )
        save_histogram(
            figure_dir / f"position_{position}_position_{target_label}.png",
            "Position [mm]",
            f"Position: position {position}",
            show,
        )
        rows.append({
            "position": position,
            "metric": "position",
            "estimator": "Position",
            "target": target,
            "event_set": "common_full_energy",
            "fit_method": "unbinned_gaussian",
            "event_count": len(valid_ids),
            "fwhm_current": position_fwhm,
            "resolution_fraction": np.nan,
            "resolution_keV": np.nan,
            "position_fwhm_mm": position_fwhm,
            "mean_current": np.nan,
            "std_current": np.nan,
        })

    if noise_values is not None:
        noise_position_values = []
        for channel_features in noise_features.values():
            if not all(np.all(np.isfinite(channel_features[channel][:2])) for channel in (0, 1)):
                continue
            ch0, ch1 = channel_features
            ratio = ch1[0] / ch0[0] if ch0[0] != 0 else np.nan
            noise_position_values.append(position_from_ratio(ratio))
        plt.figure()
        position_fwhm, _ = MakeHistgram(
            noise_position_values,
            position,
            bin_num=bin_num,
            stabilize=False,
            fit_method="unbinned_gaussian",
        )
        save_histogram(
            figure_dir / f"position_{position}_position_pulse_noise.png",
            "Position [mm]",
            f"Position: position {position}",
            show,
        )
        rows.append({
            "position": position,
            "metric": "position",
            "estimator": "Position",
            "target": "Pulse_noise",
            "event_set": "independent_noise",
            "fit_method": "unbinned_gaussian",
            "event_count": len(noise_position_values),
            "fwhm_current": position_fwhm,
            "resolution_fraction": np.nan,
            "resolution_keV": np.nan,
            "position_fwhm_mm": position_fwhm,
            "mean_current": np.nan,
            "std_current": np.nan,
        })

    result_path = output_dir / f"position_{position}_ms_vs_ms_noise.csv"
    pd.DataFrame(rows).to_csv(result_path, index=False)
    pd.DataFrame({"event_id": valid_ids}).to_csv(
        output_dir / f"position_{position}_common_event_ids.csv", index=False
    )
    return pd.DataFrame(rows)


def RunCommonAnalysis(Data_path, event_path, show=False, bin_num=None, positions=None):
    """Run the common processing for every requested absorber position."""
    data_path = Path(Data_path)
    with open(data_path / "input.json", encoding="utf-8") as file:
        para = json.load(file)
    all_positions = [int(value) for value in para["position"]]
    if positions is None:
        positions = all_positions
    else:
        positions = [int(value) for value in positions]
        unknown = sorted(set(positions) - set(all_positions))
        if unknown:
            raise ValueError(f"Positions are not in input.json: {unknown}")

    result_frames = []
    position_progress = tqdm.tqdm(
        positions,
        desc="All positions",
        position=0,
        leave=True,
        dynamic_ncols=True,
    )
    for position in position_progress:
        result_frames.append(
            CommonMSComparison(data_path, position, event_path, show, bin_num)
        )

    result_dir = data_path / "results" / "common_comparison"
    figure_dir = data_path / "figures" / "common_comparison"
    combined = pd.concat(result_frames, ignore_index=True)
    pitch = float(para["length"]) / float(para["n_abs"])
    combined["x_mm"] = (
        combined["position"].astype(float)
        - (float(para["n_abs"]) + 1.0) / 2.0
    ) * pitch
    combined_path = result_dir / f"all_positions_common_{para['E']}.csv"
    combined.to_csv(combined_path, index=False)

    # Keep the conventional files used by subScript/Show*.py.  Each target
    # receives the same common-processing data, but MS and MS+noise retain
    # only their shared FullEnergy events as indicated by event_set.
    for target in ("Pulse_ms", "Pulse_noise", "Pulse_ms_noise"):
        energy = combined[
            (combined["metric"] == "energy")
            & (combined["target"] == target)
        ]
        energy_table = (
            energy.pivot(index="x_mm", columns="estimator", values="resolution_keV")
            .reindex(columns=["Sum", "Max", "Min", "ST"])
            .sort_index()
        )
        target_result_dir = data_path / "results" / target.lower()
        target_result_dir.mkdir(parents=True, exist_ok=True)
        energy_table.to_csv(
            target_result_dir / f"ene_resos_{target}.csv",
            index_label="x_mm",
        )

        position = combined[
            (combined["metric"] == "position")
            & (combined["target"] == target)
        ].set_index("position").reindex(all_positions)
        np.savetxt(
            target_result_dir / f"fwhms_{target}.txt",
            position["position_fwhm_mm"].to_numpy(float),
        )

    target_labels = {
        "Pulse_ms": "MS",
        "Pulse_noise": "noise",
        "Pulse_ms_noise": "MS+noise",
    }
    for estimator in ("Max", "Sum"):
        figure, axis = plt.subplots(figsize=(9, 6))
        energy = combined[
            (combined["metric"] == "energy")
            & (combined["estimator"] == estimator)
        ]
        for target, label in target_labels.items():
            values = energy[energy["target"] == target].sort_values("x_mm")
            if not values.empty:
                axis.plot(values["x_mm"], values["resolution_keV"], "o-", label=label)
        axis.set_xlabel("Position [mm]")
        axis.set_ylabel("Energy Resolution [keV]")
        axis.set_title(f"{estimator}: {para['E']} keV")
        axis.grid(alpha=0.3)
        axis.legend()
        figure.tight_layout()
        figure.savefig(figure_dir / f"energy_{estimator.lower()}_all_positions_{para['E']}.png")
        plt.close(figure)

    return combined


def save_histogram(path, xlabel, title, show):
    plt.xlabel(xlabel, fontsize=15)
    plt.ylabel("Count", fontsize=15)
    plt.title(title)
    plt.tight_layout()
    plt.legend(labelspacing=0, fontsize=8, markerscale=0.5)
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(path)
    if show:
        plt.show()
    plt.close()


def _energy_column(frame):
    """Use fixed-window energy features, accepting legacy CSVs."""
    return frame["energy_height"] if "energy_height" in frame.columns else frame["height"]


def Resos(Data_path, target, show=False, bin_num=None, event_path=None):
    """Calculate position and energy resolutions from current-format feature CSVs."""
    data_path = Path(Data_path)
    if target.lower() in FULL_ENERGY_TARGETS and event_path is None:
        event_path = ask_event_path(data_path)
        if event_path is None:
            raise SystemExit(0)
    result_dir = resolution_output_dir(data_path, target)
    figure_dir = figure_output_dir(data_path, target)
    result_dir.mkdir(parents=True, exist_ok=True)
    with open(data_path / "input.json", encoding="utf-8") as f:
        para = json.load(f)
    position_from_ratio = load_ratio_calibration(data_path)
    positions = para["position"]
    colors = generate_symmetric_colors(len(positions))

    feature_pairs = []
    outlier_report = []
    quality_rows = []
    for position in positions:
        ch0, ch1, outlier_ids = load_feature_pair(
            data_path, position, target, event_path
        )
        removed = len(outlier_ids)
        if removed:
            print(f"Position {position}: removed {removed} outlier events")
        outlier_report.append({"position": position, "events": outlier_ids})
        quality_rows.append({
            "position": position,
            "outlier_or_invalid_event_count": removed,
            "resolution_event_count": len(ch0),
        })
        feature_pairs.append((ch0, ch1))

    with open(result_dir / f"reso_outliers_{target}.json", "w", encoding="utf-8") as file:
        json.dump(outlier_report, file, ensure_ascii=False, indent=2)

    feature_summary_path = result_dir / f"feature_summary_{target}.csv"
    if feature_summary_path.is_file():
        quality_summary = pd.read_csv(feature_summary_path).merge(
            pd.DataFrame(quality_rows), on="position", how="left"
        )
    else:
        quality_summary = pd.DataFrame(quality_rows)
    quality_summary.to_csv(result_dir / f"quality_summary_{target}.csv", index=False)

    position_fwhms = []
    for (ch0, ch1), position in zip(feature_pairs, positions):
        # FitRatios writes CH1 / CH0, so retain that direction here.
        ratios = ch1["height"].to_numpy(float) / ch0["height"].to_numpy(float)
        reconstructed = position_from_ratio(ratios)
        fwhm, _ = MakeHistgram(reconstructed, position, bin_num=bin_num)
        position_fwhms.append(fwhm)
    finite_position_fwhms = np.asarray(position_fwhms)[np.isfinite(position_fwhms)]
    if finite_position_fwhms.size:
        print(
            "Position resolution FWHM [mm]: "
            f"max = {np.max(finite_position_fwhms)}, "
            f"min = {np.min(finite_position_fwhms)}"
        )
    else:
        print("Position resolution FWHM [mm]: max = nan, min = nan")
    save_histogram(
        figure_dir / f"position_histgram_{target}_{para['E']}.png",
        "Position [mm]", f"{para['E']} keV", show,
    )
    np.savetxt(result_dir / f"fwhms_{target}.txt", position_fwhms)

    estimators = {
        "Sum": lambda ch0, ch1: (
            _energy_column(ch0).to_numpy(float) + _energy_column(ch1).to_numpy(float)
        ),
        "Max": lambda ch0, ch1: np.maximum(
            _energy_column(ch0).to_numpy(float), _energy_column(ch1).to_numpy(float)
        ),
        "Min": lambda ch0, ch1: np.minimum(
            _energy_column(ch0).to_numpy(float), _energy_column(ch1).to_numpy(float)
        ),
        "ST": lambda ch0, ch1: ch0["ST_Height"].to_numpy(float) + ch1["ST_Height"].to_numpy(float),
    }
    energy_resolutions = {}
    for name, estimator in estimators.items():
        resolutions = []
        for (ch0, ch1), position, color in zip(feature_pairs, positions, colors):
            _, resolution = MakeHistgram(estimator(ch0, ch1), position, color, bin_num=bin_num)
            resolutions.append(resolution * float(para["E"]))
        energy_resolutions[name] = resolutions
        finite_resolutions = np.asarray(resolutions)[np.isfinite(resolutions)]
        if finite_resolutions.size:
            print(
                f"Energy resolution {name} [keV]: "
                f"max = {np.max(finite_resolutions)}, "
                f"min = {np.min(finite_resolutions)}"
            )
        else:
            print(f"Energy resolution {name} [keV]: max = nan, min = nan")
        save_histogram(
            figure_dir / f"energy_{name.lower()}_histgram_{target}_{para['E']}.png",
            "Current [A]", f"{name}: {para['E']} keV", show,
        )

    # Match FitRatios: position is a one-based absorber-block number and the
    # centre of an odd-numbered absorber lies at x = 0.
    index = (
        np.asarray(positions, dtype=float) - (float(para["n_abs"]) + 1.0) / 2.0
    ) * float(para["length"]) / float(para["n_abs"])
    pd.DataFrame(energy_resolutions, index=index).rename_axis("x_mm").to_csv(
        result_dir / f"ene_resos_{target}.csv"
    )


def ask_data_path(initial_path=None):
    """Request a usable simulation directory when it was not passed on the CLI."""
    if initial_path is not None:
        supplied_path = Path(initial_path).expanduser()
        if (supplied_path / "input.json").is_file():
            return supplied_path
    default = str(initial_path) if initial_path else ""
    while True:
        answer = questionary.text(
            "Simulation data directory (contains input.json):", default=default
        ).ask()
        if answer is None:
            return None
        data_path = Path(answer).expanduser()
        if (data_path / "input.json").is_file():
            return data_path
        print(f"input.json was not found in: {data_path}")
        default = str(data_path)


def ask_event_path(initial_path=None):
    """Request the EventPath used by PoST_Simulation.py for external events."""
    default = str(initial_path) if initial_path else ""
    while True:
        answer = questionary.path(
            "EventPath (contains position/<id>/FullEnergyList.dat):",
            default=default,
        ).ask()
        if answer is None:
            return None
        event_path = Path(answer).expanduser()
        if event_path.is_dir():
            return event_path
        print(f"EventPath was not found or is not a directory: {event_path}")
        default = str(event_path)


def ask_mode():
    """Select the analysis stage when no mode was passed on the CLI."""
    return questionary.select(
        "Analysis mode:",
        choices=["extract", "reso", "both"],
        default="both",
    ).ask()


def ask_targets():
    """Select one target, or request every supported target interactively."""
    selection = questionary.select(
        "Pulse target:",
        choices=[*TARGET_CHOICES, ALL_TARGETS_CHOICE],
        default="Pulse_ms",
    ).ask()
    if selection is None:
        return None
    return list(TARGET_CHOICES) if selection == ALL_TARGETS_CHOICE else [selection]


def ask_bin_num():
    """Choose automatic or manually specified histogram bin count."""
    choice = questionary.select(
        "Histogram bin count:",
        choices=["automatic (optimal_bin_count)", "manual"],
    ).ask()
    if choice is None or choice.startswith("automatic"):
        return None

    while True:
        answer = questionary.text("Number of bins:", default="30").ask()
        if answer is None:
            return None
        try:
            value = int(answer)
            if value > 0:
                return value
        except ValueError:
            pass
        print("Number of bins must be a positive integer.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=(
            "Run the common TES analysis for all positions. Position ratios "
            "use adaptive 11-sample peaks; energy Max/Sum/Min use fixed "
            "noiseless-reference windows and ST uses SettlingTime."
        )
    )
    parser.add_argument(
        "data_path", nargs="?",
        help="Directory containing input.json, ratios.csv, and position folders (prompted when omitted)",
    )
    parser.add_argument(
        "--event-path",
        default=None,
        help=(
            "EventPath containing position/<id>/FullEnergyList.dat; "
            "prompted when omitted"
        ),
    )
    # Kept as a harmless compatibility flag for old command lines.  The
    # common analysis is now always the only execution path.
    parser.add_argument(
        "--common-ms-comparison", action="store_true", help=argparse.SUPPRESS
    )
    parser.add_argument(
        "--position", type=int, default=None,
        help="Analyze only this absorber block (default: all positions)",
    )
    parser.add_argument(
        "--bin-num", type=int, default=None,
        help="Display histogram bins (default: automatic; fit is always unbinned)",
    )
    parser.add_argument("--show", action="store_true", help="Display plots as well as saving them")
    arguments = parser.parse_args()
    data_path = ask_data_path(arguments.data_path)
    if data_path is None:
        raise SystemExit(0)
    event_path = arguments.event_path
    if event_path is None:
        event_path = ask_event_path(data_path)
        if event_path is None:
            raise SystemExit(0)
    if arguments.bin_num is not None and arguments.bin_num <= 0:
        parser.error("--bin-num must be positive")
    RunCommonAnalysis(
        data_path,
        event_path,
        show=arguments.show,
        bin_num=arguments.bin_num,
        positions=[arguments.position] if arguments.position is not None else None,
    )
