import glob
import os

import matplotlib.pyplot as plt
import numpy as np
import questionary
from matplotlib.widgets import RectangleSelector
from natsort import natsorted
from scipy.optimize import curve_fit
from sklearn.linear_model import RANSACRegressor

R_SH = 3.8e-3


def _display_figure(fig):
    fig.tight_layout()
    fig.canvas.draw()
    plt.show(block=True)


def _clear_figure_state():
    plt.close("all")


def CorrectJump(V_out, I_bias, idx_start, idx_stop, num_points_for_slope=10):
    V_out_range = V_out[idx_start:idx_stop]
    diff = np.abs(np.diff(V_out_range))
    jump_idx = idx_start + np.argmax(diff)

    x_points = I_bias[jump_idx - 10 : jump_idx].reshape(-1, 1)
    y_points = V_out[jump_idx - 10 : jump_idx]

    model = RANSACRegressor()
    model.fit(x_points, y_points)
    slope = model.estimator_.coef_[0]

    V_out[jump_idx + 1 :] += V_out[jump_idx] - V_out[jump_idx + 1] + slope * (
        I_bias[jump_idx + 1 :] - I_bias[jump_idx]
    )
    return V_out


def func(x, a, b):
    return a * x + b


class RangeSelector:
    def __init__(self, I_bias, V_out):
        self.I_bias = I_bias
        self.V_out = V_out
        self.selected_range = (None, None)
        self.selector = None
        self.fig = None

    def on_select(self, eclick, erelease):
        x1, x2 = sorted((eclick.xdata, erelease.xdata))
        self.selected_range = (
            np.searchsorted(self.I_bias, x1, side="left"),
            np.searchsorted(self.I_bias, x2, side="right"),
        )
        print(
            f"Selected range: {self.I_bias[self.selected_range[0]]} to "
            f"{self.I_bias[self.selected_range[1] - 1]}"
        )

    def select_range(self):
        print("Drag to select a range on the plot.")
        self.fig, ax = plt.subplots()
        ax.plot(self.I_bias, self.V_out, marker="o", c="red", linewidth=1, markersize=6)
        ax.set_title("Select Range")
        ax.set_xlabel("I_bias[uA]")
        ax.set_ylabel("V_out[V]")
        ax.grid(True)

        self.selector = RectangleSelector(
            ax,
            self.on_select,
            interactive=True,
            useblit=True,
            button=[1],
            minspanx=0,
            minspany=0,
            spancoords="data",
        )
        plt.show(block=True)
        plt.close(self.fig)
        return self.selected_range


def offset(data):
    data = np.asarray(data) - data[0]
    if np.mean(data[:5]) < 0:
        data = data * -1
    return data


def _select_path(default=None):
    prompt = f"path: [{default}]" if default else "path:"
    value = questionary.text(prompt).ask()
    return value or default


def _ensure_output_dirs(path):
    output_dir = os.path.join(path, "output")
    os.makedirs(output_dir, exist_ok=True)
    results_dir = output_dir
    return output_dir, results_dir


def _result_path(results_dir, temp):
    return os.path.join(results_dir, f"iv_{temp}.csv")


def _load_result(results_dir, temp):
    path = _result_path(results_dir, temp)
    if not os.path.exists(path):
        return None
    try:
        data = np.loadtxt(path, delimiter=",", skiprows=1)
        if data.ndim == 1:
            data = data.reshape(1, -1)
        return {
            "I_bias": data[:, 0].astype(int),
            "original": data[:, 1],
            "edited": data[:, 2],
        }
    except Exception as exc:
        print(f"Failed to load saved result for {temp}: {exc}")
        return None


def _save_result(results_dir, temp, I_bias, original, edited):
    path = _result_path(results_dir, temp)
    data = np.column_stack((I_bias, original, edited))
    np.savetxt(
        path,
        data,
        header="I_bias,Original,Edited",
        delimiter=",",
        fmt=["%d", "%.6f", "%.6f"],
    )
    return path


def _print_temperature_summary(temps, temp_files):
    print("IV temperature summary:")
    for index, (temp, files) in enumerate(zip(temps, temp_files), start=1):
        print(f"  {index}. {temp} ({len(files)} files)")


def _has_saved_results(results_dir, temps):
    for temp in temps:
        result_path = _result_path(results_dir, temp)
        if os.path.exists(result_path):
            return True
    return False


def _plot_iv_curve(ax, I_bias, v_out, label, color, marker="o", linestyle="-", markersize=6, linewidth=1):
    ax.plot(
        I_bias,
        v_out,
        marker=marker,
        linestyle=linestyle,
        linewidth=linewidth,
        markersize=markersize,
        color=color,
        label=label,
    )


def _style_iv_axes(ax, title):
    ax.set_title(title)
    ax.set_xlabel("I_bias[uA]")
    ax.set_ylabel("V_out[V]")
    ax.grid(True)
    ax.legend()


def _show_overview(temps, temp_files, results_dir):
    original_fig, original_ax = plt.subplots(figsize=(10, 6))
    saved_entries = []
    for index, (temp, files) in enumerate(zip(temps, temp_files), start=1):
        if not files:
            continue

        I_bias = []
        V_out = []
        for file_path in files:
            data = np.loadtxt(file_path)
            V_out.append(np.mean(data))
            I = os.path.splitext(os.path.basename(file_path))[0][10:-2]
            I_bias.append(int(I))

        I_bias = np.array(I_bias)
        original = offset(np.array(V_out))
        color = f"C{(index - 1) % 10}"
        _plot_iv_curve(original_ax, I_bias, original, temp, color, marker="o", linestyle="-", markersize=4, linewidth=1)

        saved = _load_result(results_dir, temp)
        if saved is not None and len(saved["edited"]) == len(I_bias):
            saved_entries.append((temp, color, saved["I_bias"], saved["edited"]))

    _style_iv_axes(original_ax, "IV overview: original")
    _display_figure(original_fig)
    plt.close(original_fig)
    _clear_figure_state()

    saved_fig = None
    if saved_entries:
        saved_fig, saved_ax = plt.subplots(figsize=(10, 6))
        for temp, color, I_bias, edited in saved_entries:
            _plot_iv_curve(saved_ax, I_bias, edited, temp, color, marker="o", linestyle="-", markersize=4, linewidth=1)
        _style_iv_axes(saved_ax, "IV overview: saved result")
        _display_figure(saved_fig)
        plt.close(saved_fig)
        _clear_figure_state()

    return original_fig, saved_fig


def _select_source_action(has_saved):
    choices = ["Edit original"]
    if has_saved:
        choices.append("Edit saved result")
    choices.append("Exit")
    return questionary.select(
        "IV analysis: choose the data source to edit or exit.",
        choices=choices,
    ).ask()


def _select_temp_action(temp, index, total, is_last):
    final_choice = "Save and Exit" if is_last else "Next temperature"
    return questionary.select(
        f"Temperature {index}/{total}: {temp}\nChoose an action",
        choices=[
            "[1] Single jump",
            "[2] Linear fit",
            "[3] Back",
            "[4] Help",
            f"[0] {final_choice}",
        ],
    ).ask()


def main(path=None):
    if path is None:
        path = _select_path()
    path = os.path.abspath(path)
    os.chdir(path)
    output_dir, results_dir = _ensure_output_dirs(path)

    temps = natsorted(glob.glob("*mK"))
    temp_files = [natsorted(glob.glob(f"{temp}/*.dat")) for temp in temps]
    _print_temperature_summary(temps, temp_files)
    has_saved = _has_saved_results(results_dir, temps)

    overview_fig = None
    saved_overview_fig = None
    if temps:
        overview_fig, saved_overview_fig = _show_overview(temps, temp_files, results_dir)

    source_choice = _select_source_action(has_saved)
    if source_choice == "Exit":
        return
    use_saved = source_choice == "Edit saved result"

    edited_records = []

    for temp_index, temp in enumerate(temps, start=1):
        files = temp_files[temp_index - 1]
        if not files:
            continue

        I_bias = []
        original_v_out = []
        for file_path in files:
            data = np.loadtxt(file_path)
            original_v_out.append(np.mean(data))
            I = os.path.splitext(os.path.basename(file_path))[0][10:-2]
            I_bias.append(int(I))

        I_bias = np.array(I_bias)
        original_v_out = offset(np.array(original_v_out))

        saved = _load_result(results_dir, temp)
        if use_saved and saved is not None and len(saved["edited"]) == len(I_bias):
            working_v_out = saved["edited"].copy()
        else:
            working_v_out = original_v_out.copy()

        previous_v_out = working_v_out.copy()
        selector = RangeSelector(I_bias, working_v_out)
        is_last = temp_index == len(temps)

        while True:
            original_fig, original_ax = plt.subplots(figsize=(10, 6))
            _plot_iv_curve(original_ax, I_bias, original_v_out, "Original", "blue")
            _plot_iv_curve(original_ax, I_bias, working_v_out, "Working", "red")
            _style_iv_axes(original_ax, f"{temp_index}/{len(temps)} {temp} original")
            _display_figure(original_fig)
            plt.close(original_fig)
            _clear_figure_state()

            saved_fig = None
            if saved is not None and len(saved["edited"]) == len(I_bias):
                saved_fig, saved_ax = plt.subplots(figsize=(10, 6))
                _plot_iv_curve(saved_ax, saved["I_bias"], saved["edited"], "Saved result", "green")
                _style_iv_axes(saved_ax, f"{temp_index}/{len(temps)} {temp} saved")
                _display_figure(saved_fig)
                plt.close(saved_fig)
                _clear_figure_state()

            choice = _select_temp_action(temp, temp_index, len(temps), is_last)
            if choice is None:
                print("No selection made. Exiting.")
                return

            if choice == "[1] Single jump":
                previous_v_out = working_v_out.copy()
                idx_start, idx_stop = selector.select_range()
                working_v_out = CorrectJump(working_v_out, I_bias, idx_start, idx_stop)
            elif choice == "[2] Linear fit":
                previous_v_out = working_v_out.copy()
                idx_start, idx_stop = selector.select_range()
                try:
                    I_range = I_bias[idx_start:idx_stop]
                    popt, cov = curve_fit(func, I_bias[:5], working_v_out[:5])
                    working_v_out[idx_start:idx_stop] = func(I_range, *popt)
                except Exception as exc:
                    print(f"Error during linear fit: {exc}")
            elif choice == "[3] Back":
                working_v_out = previous_v_out.copy()
                selector = RangeSelector(I_bias, working_v_out)
            elif choice == "[4] Help":
                print("[1] Single Jump\n局所的な段差を手で補正します。")
                print("[2] Linear fit\n選択範囲を線形に置き換えます。")
                print("[3] Back\n1つ前の状態に戻します。")
                if is_last:
                    print("[0] Save and Exit\nこの温度の補正を保存して終了します。")
                else:
                    print("[0] Next temperature\nこの温度の補正を保存して次の温度へ進みます。")
            elif choice in ("[0] Next temperature", "[0] Save and Exit"):
                edited_records.append((temp, I_bias, working_v_out.copy()))

                original_path = os.path.join(output_dir, f"iv_{temp}_original.png")
                original_fig.savefig(original_path, dpi=300, bbox_inches="tight")
                print(f"Saved {original_path}")

                if saved_fig is not None:
                    saved_path = os.path.join(output_dir, f"iv_{temp}_saved.png")
                    saved_fig.savefig(saved_path, dpi=300, bbox_inches="tight")
                    print(f"Saved {saved_path}")

                temp_result_path = _save_result(results_dir, temp, I_bias, original_v_out, working_v_out)
                print(f"Saved {temp_result_path}")

                calibration_path = os.path.join(path, "Calibration")
                os.makedirs(calibration_path, exist_ok=True)
                np.savetxt(
                    os.path.join(calibration_path, f"{temp}.dat"),
                    np.column_stack((I_bias, working_v_out)),
                    header="I_bias, V_out",
                    fmt="%f",
                )
                break

        plt.close(original_fig)
        if saved_fig is not None:
            plt.close(saved_fig)

    if overview_fig is not None:
        overview_path = os.path.join(output_dir, "iv_original_overview.png")
        overview_fig.savefig(overview_path, dpi=300, bbox_inches="tight")
        print(f"Saved {overview_path}")
        plt.close(overview_fig)
        _clear_figure_state()

    if edited_records:
        summary_fig, summary_ax = plt.subplots()
        for temp, I_bias, V_out in edited_records:
            summary_ax.plot(I_bias, V_out, label=temp, marker="o")

        summary_ax.set_xlabel("I_bias[uA]")
        summary_ax.set_ylabel("V_out[V]")
        summary_ax.set_title("I-V")
        summary_ax.grid(True)
        summary_ax.legend()
        summary_fig.savefig(os.path.join(output_dir, "iv_calib_matome.png"), dpi=300, bbox_inches="tight")
        _display_figure(summary_fig)
        plt.close(summary_fig)
        _clear_figure_state()

        ir_fig, ir_ax = plt.subplots()
        for temp, I_bias, V_out in edited_records:
            popt, cov = curve_fit(func, I_bias[:10], V_out[:10])
            eta = 1 / popt[0]
            I_tes = eta * V_out
            I_sh = I_bias - I_tes
            V_tes = I_sh * R_SH
            R_tes = np.append(0.0, V_tes[1:] / I_tes[1:])
            ir_ax.plot(I_bias, R_tes, label=temp, marker="o")

        ir_ax.set_xlabel("I_bias[uA]")
        ir_ax.set_ylabel("R[$\\Omega$]")
        ir_ax.set_title("I-R")
        ir_ax.grid(True)
        ir_ax.legend()
        ir_fig.savefig(os.path.join(output_dir, "iv_ir_calib_matome.png"), dpi=300, bbox_inches="tight")
        _display_figure(ir_fig)
        plt.close(ir_fig)
        _clear_figure_state()


if __name__ == "__main__":
    main()
