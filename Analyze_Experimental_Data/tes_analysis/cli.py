import glob
import os
import re
import shutil
import sys
from pathlib import Path

from . import analysis_utils as general
import numpy as np
import pandas as pd

from . import dispatch as analysis_dispatch
from . import operations as exp_process
from . import prompts


def _ensure_config(path):
    config_path = os.path.join(path, "PulseConfig.json")
    if not os.path.exists(config_path):
        template_path = Path(__file__).resolve().parent.parent / "PulseConfig.json"
        shutil.copy(template_path, config_path)
        print("PulseConfig.json is Copied.\nPlease set the config file and run again.")
        sys.exit(0)
    return general.LoadJson(config_path)


def _get_channels(path):
    folders = glob.glob(os.path.join(path, "CH*_pulse"))
    return [re.search(r"CH(.*)_pulse", os.path.basename(folder)).group(1) for folder in folders]


def _select_mode():
    return prompts.select_mode()


def _load_pulse_output(path, channel):
    """output.csvを読み、過去の最適フィルタ結果があればPeakOptも軸候補に加える。"""
    folder = os.path.join(path, f"CH{channel}_pulse")
    df = pd.read_csv(os.path.join(folder, "output.csv"))

    opt_path = os.path.join(folder, "output_optimalfilter.csv")
    if os.path.exists(opt_path) and "PeakOpt" not in df.columns:
        opt = pd.read_csv(opt_path)
        if "PeakOpt" in opt.columns:
            # 前回の実行で絞り込まれたキーにしか値が無いので、残りはNaNになる。
            df = df.merge(opt[["key", "PeakOpt"]], on="key", how="left")
            print(f"CH{channel}: 前回のPeakOptを軸候補に追加しました（{opt_path}）")
    return df


def _axis_choices(*dfs):
    """散布図の軸に使える数値列を、全DataFrameに共通するものだけ返す。"""
    choices = None
    for df in dfs:
        numeric = [
            col for col in df.columns
            if col != "key" and pd.api.types.is_numeric_dtype(df[col])
        ]
        if choices is None:
            choices = numeric
        else:
            choices = [col for col in choices if col in numeric]
    return choices or list(prompts.KEY_CHOICES)


def _keep_selected(df, selected_keys):
    if selected_keys is None:
        return df
    return df[df["key"].isin(selected_keys)].reset_index(drop=True)


def _refine_done(selected_keys):
    """絞り込みを続けるか判定する。Trueなら選択を確定して抜ける。"""
    count = len(selected_keys)
    print(f"選択済み: {count} 件")
    return prompts.select_refine_action(count) != "Select again"


def _pick_single(df, axes, selected_keys=None):
    x_key = prompts.select_x_key(axes)
    y_key = prompts.select_y_key(axes)
    # selected_keysがあれば、その中だけを対象にしてさらに絞り込む。
    return general.SelectIDFrom1DF(_keep_selected(df, selected_keys), x_key, y_key)


def _pick_two(df1, df2, axes, selected_keys=None):
    key = prompts.select_key(axes)
    picked = general.SelectIDFrom2DF(
        _keep_selected(df1, selected_keys),
        _keep_selected(df2, selected_keys),
        key,
    )
    return picked, key


def _choose_two_channels(chs):
    if len(chs) < 2:
        print("2チャンネル以上のデータが必要です。")
        return None
    if len(chs) == 2:
        return list(chs)

    channels = prompts.select_two_channels(chs)
    if not channels:
        return None
    if len(channels) != 2:
        print("2チャンネルを選択してください。")
        return None
    return channels


def _select_single_channel_keys(path, chs):
    channel = prompts.select_channel(chs)
    df = _load_pulse_output(path, channel)
    return channel, _pick_single(df, _axis_choices(df))


def _select_two_channel_keys(path, chs):
    channels = _choose_two_channels(chs)
    if channels is None:
        return None, None, None

    df1 = _load_pulse_output(path, channels[0])
    df2 = _load_pulse_output(path, channels[1])
    selected_keys, key = _pick_two(df1, df2, _axis_choices(df1, df2))
    return channels, key, selected_keys


def _save_selected_keys(path, selected_keys):
    output_dir = os.path.join(path, "output")
    os.makedirs(output_dir, exist_ok=True)
    output_path = os.path.join(output_dir, "pulse_SelectedKeys.txt")
    np.savetxt(output_path, selected_keys, fmt="%d")
    print(f"Saved selected keys to {output_path}")


def _load_noise_models(path, channels):
    noises = {}
    for ch in channels:
        noise_path = exp_process.NoiseModelPath(path, ch)
        if not os.path.exists(noise_path):
            raise FileNotFoundError(
                f"Noise model not found: {noise_path}. "
                "Run Noise Analysis first."
            )
        noises[ch] = general.LoadTxt(noise_path)
    return noises


def _input_eta():
    text = prompts.input_eta()
    try:
        return float(text)
    except (TypeError, ValueError):
        print("etaは数値で入力してください。")
        return None


def _run_temp_and_optimal(config, path):
    chs = _get_channels(path)
    mode = _select_mode()

    if mode == "Single Channel":
        channels = [prompts.select_channel(chs)]
    elif mode == "Two Channels":
        channels = _choose_two_channels(chs)
        if channels is None:
            return
    else:
        return

    noises = _load_noise_models(path, channels)
    # 選択のたびに最適フィルタを掛け直すので、etaは最初に一度だけ聞く。
    eta = _input_eta()
    if eta is None:
        return

    selected_keys = None
    while True:
        # 直前のラウンドで更新されたPeakOptを軸候補に含めるため、毎回読み直す。
        dfs = [_load_pulse_output(path, ch) for ch in channels]
        axes = _axis_choices(*dfs)

        if len(channels) == 1:
            picked = _pick_single(dfs[0], axes, selected_keys)
        else:
            picked, _key = _pick_two(dfs[0], dfs[1], axes, selected_keys)

        if len(picked) > 0:
            selected_keys = picked
        else:
            print("選択されたデータがありません。直前の選択を維持します。")

        if selected_keys is None or len(selected_keys) == 0:
            print("選択されたデータがありません。処理を中止します。")
            return

        _save_selected_keys(path, selected_keys)

        # 選択するたびに最適フィルタを作り直し、PeakOptを更新する。
        for ch in channels:
            exp_process.OptimalFilter(
                config, path, noises[ch], ch, selected_keys,
                eta_uA_per_V=eta,
            )

        if _refine_done(selected_keys):
            break

    # TempCalibは全CH*_pulseを走査するので、全チャンネルの
    # output_optimalfilter.csvが出そろってから一度だけ実行する。
    exp_process.TempCalib(path, selected_keys)


def _run_select_from_scatter(path):
    chs = _get_channels(path)
    mode = _select_mode()

    if mode == "Single Channel":
        _channel, selected_keys = _select_single_channel_keys(path, chs)
        output_dir = os.path.join(path, "output")
        os.makedirs(output_dir, exist_ok=True)
        np.savetxt(os.path.join(output_dir, "pulse_SelectedKeys.txt"), selected_keys, fmt="%d")
    elif mode == "Two Channels":
        channels, _key, selected_keys = _select_two_channel_keys(path, chs)
        if channels is None:
            return
        output_dir = os.path.join(path, "output")
        os.makedirs(output_dir, exist_ok=True)
        np.savetxt(os.path.join(output_dir, "pulse_SelectedKeys.txt"), selected_keys, fmt="%d")

    print(f"Saved selected keys to {os.path.join(path, 'output', 'pulse_SelectedKeys.txt')}")


def _run_compare_estimators(config, path):
    """PeakOpt系カラムの分解能をまとめて比較する（項目7・8）。"""
    chs = _get_channels(path)
    channel = prompts.select_channel(chs)
    folder = os.path.join(path, f"CH{channel}_pulse")

    candidates = [
        name for name in ("output_tempcalib.csv", "output_optimalfilter.csv")
        if os.path.exists(os.path.join(folder, name))
    ]
    csvs = [os.path.basename(p) for p in glob.glob(os.path.join(folder, "*.csv"))]
    files = candidates + [name for name in csvs if name not in candidates]
    if not files:
        print(f"{folder} に比較できるCSVがありません。")
        return

    file = prompts.select_csv_file(files)
    csv_path = os.path.join(folder, file)

    exp_process.CompareEstimators(csv_path)
    exp_process.CompareChain(config, path, channel, LoadPath=file)


def _run_view_pulse(path):
    chs = _get_channels(path)
    channel = prompts.select_channel(chs)
    key = prompts.input_integer("Input Key (integer):")
    try:
        key = int(key)
    except ValueError:
        print("Keyは整数で入力してください。")
        return

    exp_process.ViewPulse(path, channel, key)


def _run_hist(path):
    chs = _get_channels(path)
    channel = prompts.select_channel(chs)

    csvs = glob.glob(os.path.join(path, f"CH{channel}_pulse", "*.csv"))
    files = [os.path.basename(csv_path) for csv_path in csvs]
    file = prompts.select_csv_file(files)

    csv_path = os.path.join(path, f"CH{channel}_pulse", file)
    csv = pd.read_csv(csv_path)
    keys = list(csv.columns)
    key = prompts.select_csv_column(keys)

    bin_choose = prompts.select_bin_option()
    if bin_choose == "Manual":
        bin_num = prompts.input_integer("Input bin number (integer):")
        try:
            bin_num = int(bin_num)
        except ValueError:
            print("Bin numberは整数で入力してください。")
            return
        exp_process.Hist(csv_path, key, binNum=bin_num)
        return

    exp_process.Hist(csv_path, key)


def pulse_main(path):
    config = _ensure_config(path)

    choice = prompts.select_analysis_type()
    if choice == "Pulse Analysis":
        exp_process.PulseAnalysis(config, path)
    elif choice == "Noise Analysis":
        exp_process.NoiseAnalysis(config, path)
    elif choice == "Temp and Optimal":
        _run_temp_and_optimal(config, path)
    elif choice == "Compare Estimators":
        _run_compare_estimators(config, path)
    elif choice == "Scatter2D":
        exp_process.Scatter2D(path)
    elif choice == "Select from Scatter":
        _run_select_from_scatter(path)
    elif choice == "ViewPulse":
        _run_view_pulse(path)
    elif choice == "Hist":
        _run_hist(path)
    elif choice == "Exit":
        sys.exit(0)


def _run_pulse_loop(path):
    while True:
        pulse_main(path)


def run():
    path = general.InputPath()
    analysis_dispatch.dispatch(path, _run_pulse_loop)


if __name__ == "__main__":
    run()
