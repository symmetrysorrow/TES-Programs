"""Interactive viewer for current post_all HDF5 pulse files.

Examples
--------
python observe.py H:\\hata2025\\New_test

The reader is streaming: selecting one event or creating a histogram does
not load the multi-gigabyte JSON document into memory.
"""

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import questionary
from matplotlib.widgets import RectangleSelector
from tqdm import tqdm

from lib.pulse_hdf5 import iter_pulse_items, read_pulse


MODE_FILES = {
    "normal": "pulses.h5",
    "noise": "pulse_noise.h5",
    "ms": "pulse_MS.h5",
    "noise_ms": "pulse_MS_noise.h5",
}

FEATURE_FILES = {
    "normal": ("pulses_output_TES0.csv", "pulses_output_TES1.csv"),
    "noise": ("pulse_noise_output_TES0.csv", "pulse_noise_output_TES1.csv"),
    "ms": ("pulse_MS_output_TES0.csv", "pulse_MS_output_TES1.csv"),
    "noise_ms": ("pulse_MS_noise_output_TES0.csv", "pulse_MS_noise_output_TES1.csv"),
}


def load_input(data_path):
    with open(Path(data_path) / "input.json", encoding="utf-8") as file:
        return json.load(file)


def pulse_path(data_path, mode, position):
    data_path = Path(data_path)
    if mode == "normal":
        return data_path / MODE_FILES[mode]
    return data_path / str(position) / MODE_FILES[mode]


def event_id_for(mode, position, requested_id):
    # pulses.h5 has absorber positions as its event IDs, not PHITS event IDs.
    return str(position) if mode == "normal" else str(requested_id)


def find_pulse(path, event_id):
    """Return one pulse directly from the HDF5 event axis."""
    return read_pulse(path, event_id)


def plot_waveform(pulse, rate, title):
    try:
        ch0 = np.asarray(pulse["ch0"], dtype=float)
        ch1 = np.asarray(pulse["ch1"], dtype=float)
    except (KeyError, TypeError) as error:
        raise ValueError("Pulse must contain numeric ch0 and ch1 arrays") from error
    if len(ch0) != len(ch1):
        raise ValueError(f"Channel length mismatch: CH0={len(ch0)}, CH1={len(ch1)}")

    time_ms = np.arange(len(ch0)) / float(rate) * 1e3
    figure, axis = plt.subplots(figsize=(11, 5))
    axis.plot(time_ms, ch0, label="CH0", linewidth=1.0)
    axis.plot(time_ms, ch1, label="CH1", linewidth=1.0)
    axis.set(title=title, xlabel="Time [ms]", ylabel="Current [A]")
    axis.grid(alpha=0.3)
    axis.legend()
    figure.tight_layout()
    plt.show()


METRICS = {
    "CH0 peak": lambda ch0, ch1: np.max(ch0),
    "CH1 peak": lambda ch0, ch1: np.max(ch1),
    "CH0 + CH1 peak": lambda ch0, ch1: np.max(ch0) + np.max(ch1),
    "CH1 / CH0 peak": lambda ch0, ch1: np.max(ch1) / np.max(ch0),
    "CH0 peak index": lambda ch0, ch1: int(np.argmax(ch0)),
}


def histogram_values(path, metric_name):
    """Extract one scalar per event while retaining only the scalar array."""
    metric = METRICS[metric_name]
    values = []
    for event_id, pulse in tqdm(iter_pulse_items(path), desc=f"Reading {path.name}", unit="event"):
        try:
            ch0 = np.asarray(pulse["ch0"], dtype=float)
            ch1 = np.asarray(pulse["ch1"], dtype=float)
            value = metric(ch0, ch1)
        except (KeyError, TypeError, ValueError, FloatingPointError) as error:
            print(f"Skipping invalid event {event_id}: {error}")
            continue
        if np.isfinite(value):
            values.append(value)
    return np.asarray(values, dtype=float)


def plot_histogram(values, metric_name, title):
    if len(values) == 0:
        raise ValueError("No finite values were available for the histogram")
    figure, axis = plt.subplots(figsize=(9, 5))
    axis.hist(values, bins="auto", color="tab:blue", alpha=0.85)
    axis.set(title=title, xlabel=metric_name, ylabel="Count")
    axis.grid(alpha=0.3)
    figure.tight_layout()
    plt.show()


def plot_feature_histogram(data_path, mode, position):
    """Plot a feature CSV histogram using all or selected event IDs."""
    folder = Path(data_path) / str(position)
    file_names = FEATURE_FILES[mode]
    channel = questionary.select("Channel:", choices=["CH0", "CH1", "Sum"]).ask()
    target = questionary.select(
        "Feature:", choices=["height", "rise", "ST_Height"]
    ).ask()
    source = questionary.select(
        "Events:", choices=["all", "selected_ids.txt"]
    ).ask()
    if None in (channel, target, source):
        return

    if channel == "Sum":
        frame0 = pd.read_csv(folder / file_names[0], index_col="id")
        frame1 = pd.read_csv(folder / file_names[1], index_col="id")
        common_ids = frame0.index.intersection(frame1.index, sort=False)
        if target not in frame0 or target not in frame1:
            raise ValueError(f"Feature column was not found: {target}")
        frame = pd.DataFrame(
            {target: frame0.loc[common_ids, target] + frame1.loc[common_ids, target]},
            index=common_ids,
        )
    else:
        csv_path = folder / file_names[0 if channel == "CH0" else 1]
        if not csv_path.is_file():
            raise FileNotFoundError(f"Feature CSV was not found: {csv_path}")
        frame = pd.read_csv(csv_path, index_col="id")
        if target not in frame:
            raise ValueError(f"Feature column was not found: {target}")

    if source == "selected_ids.txt":
        ids_path = folder / "selected_ids.txt"
        if not ids_path.is_file():
            raise FileNotFoundError(f"Selected ID file was not found: {ids_path}")
        selected_ids = ids_path.read_text(encoding="utf-8").splitlines()
        frame = frame.loc[frame.index.astype(str).isin(selected_ids)]

    values = frame[target].to_numpy(float)
    values = values[np.isfinite(values)]
    if len(values) == 0:
        raise ValueError("No finite feature values are available")
    figure, axis = plt.subplots(figsize=(9, 5))
    axis.hist(values, bins="auto", color="tab:blue", alpha=0.85)
    axis.set(
        title=f"{channel} {target}, {source}, position {position}",
        xlabel=target, ylabel="Count",
    )
    axis.grid(alpha=0.3)
    figure.tight_layout()
    plt.show()


def plot_scatter(data_path, mode, position):
    """Show a selectable scatter plot from current feature CSV files."""
    folder = Path(data_path) / str(position)
    file0, file1 = (folder / name for name in FEATURE_FILES[mode])
    if not file0.is_file() or not file1.is_file():
        raise FileNotFoundError(
            f"Feature CSVs were not found in {folder}. Run reso.py extract first."
        )

    df0 = pd.read_csv(file0, index_col="id")
    df1 = pd.read_csv(file1, index_col="id")
    common_ids = df0.index.intersection(df1.index, sort=False)
    columns = ["height", "rise", "ST_Height"]
    missing = [column for column in columns if column not in df0 or column not in df1]
    if missing:
        raise ValueError(f"Feature CSVs are missing columns: {', '.join(missing)}")

    x_channel = questionary.select("X axis channel:", choices=["CH0", "CH1"]).ask()
    x_target = questionary.select("X axis value:", choices=columns).ask()
    y_channel = questionary.select("Y axis channel:", choices=["CH0", "CH1"]).ask()
    y_target = questionary.select("Y axis value:", choices=columns).ask()
    if None in (x_channel, x_target, y_channel, y_target):
        return

    x_df = df0 if x_channel == "CH0" else df1
    y_df = df0 if y_channel == "CH0" else df1
    x_data = x_df.loc[common_ids, x_target].to_numpy(float)
    y_data = y_df.loc[common_ids, y_target].to_numpy(float)
    ids = common_ids.to_numpy()
    valid = np.isfinite(x_data) & np.isfinite(y_data)
    x_data, y_data, ids = x_data[valid], y_data[valid], ids[valid]
    if len(ids) == 0:
        raise ValueError("No finite feature pairs are available for the scatter plot")

    figure, axis = plt.subplots(figsize=(9, 7))
    axis.scatter(x_data, y_data, s=8, alpha=0.6)
    selected_ids = set()

    def onselect(eclick, erelease):
        if eclick.xdata is None or erelease.xdata is None or eclick.ydata is None or erelease.ydata is None:
            return
        x_min, x_max = sorted((eclick.xdata, erelease.xdata))
        y_min, y_max = sorted((eclick.ydata, erelease.ydata))
        selected = (x_min <= x_data) & (x_data <= x_max) & (y_min <= y_data) & (y_data <= y_max)
        selected_ids.update(str(event_id) for event_id in ids[selected])
        output = folder / "selected_ids.txt"
        output.write_text("".join(f"{event_id}\n" for event_id in sorted(selected_ids)), encoding="utf-8")
        print(f"Saved {len(selected_ids)} IDs to {output}")
        axis.scatter(x_data[selected], y_data[selected], color="red", s=12)
        figure.canvas.draw_idle()

    selector = RectangleSelector(
        axis, onselect, useblit=True, button=[1], interactive=True
    )
    figure._scatter_selector = selector
    axis.set(title=f"2D Scatter: {mode}, position {position}", xlabel=f"{x_channel} {x_target}", ylabel=f"{y_channel} {y_target}")
    axis.grid(alpha=0.3)
    figure.tight_layout()
    plt.show()


def ask_position(positions):
    default = str(positions[0]) if positions else "1"
    answer = questionary.text("Position (absorber block number):", default=default).ask()
    if answer is None:
        return None
    try:
        position = int(answer)
    except ValueError:
        print("Position must be an integer.")
        return None
    if positions and position not in positions:
        print(f"Position {position} is not in input.json: {positions}")
        return None
    return position


def choose_source(data_path, event_positions, normal_positions):
    mode = questionary.select("Pulse type:", choices=list(MODE_FILES)).ask()
    if mode is None:
        return None, None
    position = ask_position(normal_positions if mode == "normal" else event_positions)
    if position is None:
        return None, None
    path = pulse_path(data_path, mode, position)
    if not path.is_file():
        print(f"File was not found: {path}")
        return None, None
    return mode, position


def interactive_main(data_path):
    parameters = load_input(data_path)
    event_positions = [int(position) for position in parameters.get("position", [])]
    normal_positions = list(range(1, int(parameters["n_abs"]) + 1))
    rate = float(parameters["rate"])

    while True:
        action = questionary.select(
            "Observation:",
            choices=[
                "View one waveform", "Show pulse histogram",
                "Show feature histogram", "Show 2D scatter", "Quit",
            ],
        ).ask()
        if action in (None, "Quit"):
            return

        mode, position = choose_source(data_path, event_positions, normal_positions)
        if mode is None:
            continue
        path = pulse_path(data_path, mode, position)

        if action == "View one waveform":
            requested_id = questionary.text("Event ID:", default="847892").ask()
            if requested_id is None:
                continue
            event_id = event_id_for(mode, position, requested_id)
            try:
                pulse = find_pulse(path, event_id)
                plot_waveform(pulse, rate, f"{mode}, position {position}, event {event_id}")
            except (KeyError, ValueError) as error:
                print(error)
        elif action == "Show pulse histogram":
            metric_name = questionary.select("Histogram value:", choices=list(METRICS)).ask()
            if metric_name is None:
                continue
            values = histogram_values(path, metric_name)
            try:
                plot_histogram(values, metric_name, f"{mode}, position {position}")
            except ValueError as error:
                print(error)
        elif action == "Show feature histogram":
            try:
                plot_feature_histogram(data_path, mode, position)
            except (FileNotFoundError, ValueError) as error:
                print(error)
        else:
            try:
                plot_scatter(data_path, mode, position)
            except (FileNotFoundError, ValueError) as error:
                print(error)


def ask_data_path(initial_path=None):
    """Request a usable simulation directory when it was not passed on the CLI."""
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


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Interactively inspect post_all HDF5 pulse files.")
    parser.add_argument(
        "data_path", nargs="?",
        help="Directory containing input.json and pulse JSON files (prompted when omitted)",
    )
    arguments = parser.parse_args()
    data_path = ask_data_path(arguments.data_path)
    if data_path is not None:
        interactive_main(data_path)
