"""Analysis operations for pulse, noise, calibration, and RT datasets.

The interactive menu lives in :mod:`pulse`; this module contains the actual
data-processing operations so they can also be called from scripts.
"""

import glob
import os
import re
import shutil

import natsort
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import questionary
import tqdm
from scipy.optimize import curve_fit
from contextlib import contextmanager

from . import analysis_utils as general
from .noise_utils import one_sided_asd_from_power, voltage_asd_to_pA


@contextmanager
def _cd(path):
    prev = os.getcwd()
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(prev)


def _load_dat(path):
    try:
        return np.loadtxt(path, comments="#", skiprows=6)
    except Exception:
        return np.loadtxt(path, comments="#")


def _linear(x, a, b):
    return a * x + b


def _offset(data):
    data = np.asarray(data) - data[0]
    if np.mean(data[:10]) < 0:
        data = data * -1
    return data


def _extract_int(pattern, text):
    match = re.search(pattern, text)
    if match is None:
        raise ValueError(f"Could not parse value from {text}")
    return int(match.group(1))


# ---------------------------------------------------------------------------
# Pulse and noise analysis
# ---------------------------------------------------------------------------

def PulseAnalysis(config: dict, path: str):
    folders=glob.glob(f"{path}/CH*_pulse")

    for folder in tqdm.tqdm(folders, desc="Pulse analysis (channels)"):
        if os.path.exists(f"{folder}/output.csv"):
            Skip=questionary.confirm(f"output.csv already exists in {folder}. Do you want to skip Pulse Analysis for this folder?").ask()
            if Skip:
                print(f"Skipped Pulse Analysis for {folder}.")
                continue

        pulse_pathes = glob.glob(f"{folder}/rawdata/CH*.dat")
        pulse_pathes = natsort.natsorted(pulse_pathes)
        results = []

        for pulse_path in tqdm.tqdm(
            pulse_pathes, desc=f"Pulses: {os.path.basename(folder)}"
        ):
            # ファイル名からキーを抽出
            filename = os.path.basename(pulse_path)
            key = os.path.splitext(filename)[0].split("_")[-1]
            
            pulse = general.LoadBin(pulse_path)
            if len(pulse) != config["Readout"]["Sample"]:
                continue

            result = general.AnalyzePulse(pulse, config,key)
            if result is not None:
                results.append(result)
            

        df = pd.DataFrame(results)
        if "key" in df.columns:
            df = df.sort_values("key").reset_index(drop=True)

        # 保存
        output_path = f"{folder}/output.csv"
        df.to_csv(output_path, index=False)
        print(f"Saved results to {output_path}")

    dfs = {}
    for folder in tqdm.tqdm(folders, desc="Synchronizing pulse keys"):
        output_path = f"{folder}/output.csv"
        if os.path.exists(output_path):
            df = pd.read_csv(output_path)
            if "key" in df.columns:
                dfs[folder] = df

    all_keys = [set(df["key"]) for df in dfs.values()]
    common_keys = set.intersection(*all_keys) if all_keys else set()

    for folder, df in dfs.items():
        df = df[df["key"].isin(common_keys)].reset_index(drop=True)
        output_path = f"{folder}/output.csv"
        df.to_csv(output_path, index=False)

def _is_legacy_filter_method(filter_method):
    return filter_method == general.LEGACY_FFT_METHOD


# 方式ごとに出力名を分けておかないと、比較のために方式を切り替えたときに
# 直前のテンプレートを上書きしてしまう。
_METHOD_SUFFIX = {
    general.LEGACY_FFT_METHOD: "_old",
    general.CURRENT_METHOD: "_new",
    general.CURRENT_NO_TEMPLATE_BESSEL_METHOD: "_nobessel",
    general.PSD_OPTIMAL_METHOD: "_psd",
}


def _remove_dc(data):
    """Remove the per-record DC component before filtering and FFT."""
    data = np.asarray(data)
    return data - np.mean(data)


def _legacy_output_name(config):
    """Return the filename used for the legacy numerical result."""
    return "modelnoise_old.txt"


def NoiseModelPath(path, channel, filter_method="Current (rfft/irfft + Bessel)", config=None):
    """Return the modelnoise path corresponding to the selected FFT policy."""
    noise_folder = os.path.join(path, f"CH{channel}_noise")
    if _is_legacy_filter_method(filter_method):
        return os.path.join(noise_folder, _legacy_output_name(config))
    return os.path.join(noise_folder, "modelnoise.txt")


def OptimalFilterPath(path, channel, filter_method="Current (rfft/irfft + Bessel)", name="opt_template", ext="txt"):
    """Return the optimal-filter output path for the selected FFT policy."""
    # getparaが出力した既存ファイルを上書きしないよう、必ずサフィックスを付ける。
    pulse_folder = os.path.join(path, f"CH{channel}_pulse")
    suffix = _METHOD_SUFFIX.get(filter_method, "_new")
    return os.path.join(pulse_folder, f"{name}{suffix}.{ext}")


def _legacy_noise_model(config, noise_paths):
    """Reproduce the numerical model used by the legacy noise_main.py."""
    sample = config["Readout"]["Sample"]
    rate = config["Readout"]["Rate"]
    cutoff = config["Analysis"]["CutoffFrequency"]
    presample = config["Readout"].get("PreSample", 0)
    threshold = config["Analysis"].get(
        "NoiseThreshold",
        config.get("Config", {}).get("threshold", np.inf),
    )

    model = np.zeros(sample)
    for noise_path in tqdm.tqdm(
        natsort.natsorted(noise_paths), desc="Legacy noise model", leave=False
    ):
        try:
            data = general.LoadBin(noise_path)
            if len(data) != sample:
                continue

            base = np.mean(data[:presample]) if presample else 0
            data_baseline = data - base
            data = _remove_dc(data)
            if cutoff > 0:
                data = general.Bessel(data, rate, cutoff)

            # This follows the legacy acceptance test.  The old baseline
            # condition was contradictory (<= -3 and >= 3), so it is kept
            # exactly for numerical compatibility.
            peak = np.max(data_baseline)
            if (base <= -3 and base >= 3) or peak >= threshold:
                continue

            model += np.abs(np.fft.fft(data))
        except (OSError, TypeError, ValueError):
            continue

    if not noise_paths:
        return None

    model /= len(noise_paths)
    df = rate / sample
    amp_dens = np.sqrt(model**2 / df)
    return amp_dens[: sample // 2 + 1]


def NoiseAnalysis(
    config: dict,
    path: str,
    FilterMethod="Current (rfft/irfft + Bessel)",
    eta_uA_per_V=None,
):
    
    # 設定値の取得
    sample = config["Readout"]["Sample"] # データ点数 (N)
    rate = config["Readout"]["Rate"]     # サンプリングレート (Fs)
    cutoff = config["Analysis"]["CutoffFrequency"] # ローパスフィルタのカットオフ周波数

    if eta_uA_per_V is None:
        eta_uA_per_V = float(input("eta [uA/V]:"))
    else:
        eta_uA_per_V = float(eta_uA_per_V)

    noise_threshold=0.04
    
    # フォルダごとの処理
    for folder in tqdm.tqdm(
        glob.glob(f"{path}/CH*_noise"), desc="Noise analysis (channels)"
    ):
        noise_pathes = glob.glob(f"{folder}/rawdata/CH*.dat")

        if _is_legacy_filter_method(FilterMethod):
            amp_dens = _legacy_noise_model(config, noise_pathes)
            if amp_dens is None:
                print(f"No noise data found in {folder}; skipping.")
                continue

            noise_model_path = NoiseModelPath(
                path,
                os.path.basename(folder).removeprefix("CH").removesuffix("_noise"),
                FilterMethod,
                config,
            )
            amp_dens_pA = voltage_asd_to_pA(amp_dens, eta_uA_per_V)
            np.savetxt(
                noise_model_path,
                amp_dens_pA,
                header=(
                    "one-sided amplitude spectral DENSITY (ASD) [pA/sqrt(Hz)];"
                    " PSD = this**2 ; legacy amplitude-averaged model"
                ),
            )
            plt.plot(
                general.GetFreq(rate, sample)[: sample // 2 + 1],
                amp_dens_pA,
                linestyle="-",
                linewidth=0.7,
            )
            plt.loglog()
            plt.xlabel("Frequency[Hz]")
            plt.ylabel("Intensity[pA/Hz$^{1/2}$]")
            plt.grid()
            plt.savefig(os.path.splitext(noise_model_path)[0] + ".png")
            plt.show()
            continue
        
        power_model = np.zeros(sample // 2 + 1)
        count = 0
        original_noise_list = []

        # メディアンフィルタのカーネルサイズ (窓サイズ) を定義
        # 奇数に設定し、ノイズの幅に応じて調整します。3, 5, 7 などが一般的です。
        # ここでは例としてカーネルサイズ3を使用します。
        median_kernel_size = 3

        # 窓関数なしでFFTすると、矩形窓のサイドローブを通じて低周波成分が
        # 高周波ビンへ漏れ込み(スペクトル漏れ)、実在しない高周波ノイズ床が
        # 見かけ上生じる。Hann窓でこれを抑える。
        # ENBW補正(sqrt(mean(window**2)))を使うのは、対象が正弦波ではなく
        # 広帯域ランダムノイズのパワーだから(コヒーレントゲイン補正ではない)。
        fft_window = np.hanning(sample)
        window_power_gain = np.sqrt(np.mean(fft_window ** 2))

        filtered=0
        sample_unmatch=0
        
        for noise_path in tqdm.tqdm(
            noise_pathes, desc=f"Noise: {os.path.basename(folder)}"
        ):
            noise = general.LoadBin(noise_path)
            if len(noise) != sample:
                sample_unmatch+=1
                continue

            diff=np.max(noise)-np.min(noise)
            if diff>noise_threshold:
                filtered+=1
                continue
            
            # 1. スパイクノイズ除去（メディアンフィルタ）を追加
            # scipy.signal.medfilt を使用するには、事前に import scipy.signal が必要です。
            #noise = scipy.signal.medfilt(noise, kernel_size=median_kernel_size)
            
            # 2. 既存の処理を続行
            # 各レコードの大きなDCオフセットをHann窓へ入れると、窓の有限長
            # スペクトルが高周波側まで漏れ込み、特にデジタルBesselがほぼ0に
            # なるNyquist近傍で偽のノイズ床として支配的になる。平均値だけを
            # 引けばACノイズ成分は保ったまま、この解析アーティファクトを除ける。
            noise = _remove_dc(noise)
            noise = general.Bessel(noise, rate, cutoff)
            diff=np.max(noise)-np.min(noise)
            if diff>noise_threshold:
                filtered+=1
                continue
            original_noise_list.append(noise)

            noise_fft = np.fft.rfft(noise * fft_window)
            power_model += np.abs(noise_fft) ** 2
            count += 1

        print(f"count;{count}")

        print(f"Filtered {filtered} and skipped {sample_unmatch} from {len(noise_pathes)}.")

        if count == 0:
            print(f"No usable noise data found in {folder}; skipping.")
            continue

        mean_power = power_model / count
        # Keep the historical filename for downstream users, but make its
        # contents explicit: RMS (power-averaged) Hann-window FFT magnitude
        # before ASD and eta normalization.
        np.savetxt(
            f"{folder}/noise_fft_Amplitude.txt",
            np.sqrt(mean_power),
            header="RMS magnitude sqrt(mean(abs(rfft(noise * Hann))**2))",
        )

        amp_dens = one_sided_asd_from_power(
            mean_power,
            sample,
            rate,
            window_power_gain,
        )
        amp_dens = voltage_asd_to_pA(amp_dens, eta_uA_per_V)

        noise_model_path = NoiseModelPath(
            path,
            os.path.basename(folder).removeprefix("CH").removesuffix("_noise"),
            FilterMethod,
            config,
        )
        os.makedirs(os.path.dirname(noise_model_path), exist_ok=True)
        # 保存しているのは ASD（振幅密度）であってPSDでもFFT振幅でもない。
        # 最適フィルタで S*/PSD を組むときは PSD = ASD**2 を使うこと。
        np.savetxt(
            noise_model_path,
            amp_dens,
            header=(
                "one-sided amplitude spectral DENSITY (ASD) [pA/sqrt(Hz)];"
                " PSD = this**2 ; power-averaged, Hann window, ENBW corrected"
            ),
        )

        fq=general.GetFreq(rate,sample)

        # スペクトルをグラフ化
        plt.plot(fq[: int(sample / 2) + 1], amp_dens, linestyle="-", linewidth=0.7)
        plt.loglog()
        plt.xlabel("Frequency[Hz]")
        plt.ylabel("Intensity[pA/Hz$^{1/2}$]")
        plt.grid()
        plt.savefig(os.path.splitext(noise_model_path)[0] + ".png")
        plt.show()
        


# ---------------------------------------------------------------------------
# Calibration and optimal filtering
# ---------------------------------------------------------------------------

# 波高推定量として比較したいカラム（存在するものだけ使う）。
ESTIMATOR_COLUMNS = (
    "Peak",
    "PeakOptLegacy",
    "PeakOpt",
    "PeakOptPSD",
)

# TempCalib 後のカラム名（推定量 -> 補正後）。
CALIBRATED_SUFFIX = "Temp"


def _estimator_columns(df):
    return [col for col in ESTIMATOR_COLUMNS if col in df.columns]


def _calibrated_name(column):
    # PeakOpt -> PeakOptTemp（既存の後段処理が使う名前を維持する）。
    return f"{column}{CALIBRATED_SUFFIX}"


def TempCalib(
    path: str,
    SelectedKeys,
    SavePath="output_tempcalib.csv",
    LoadPath="output_optimalfilter.csv",
    Columns=None,
    plot=True,
):
    """全推定量カラムについてベースライン補正を掛け、分解能を比較表示する。

    ``PeakOpt`` -> ``PeakOptTemp`` は従来どおり。加えて ``Peak`` や
    ``PeakOptPSD`` も同じ手続きで補正するので、どの段階で Getpara と差が
    出るのかを補正前後の両方で追える。
    """
    results = {}
    for folder in tqdm.tqdm(glob.glob(f"{path}/CH*_pulse")):
        if not os.path.exists(f"{folder}/{LoadPath}"):
            print(f"output.csv not found in {folder}/{LoadPath}, skipping TempCalib.")
            continue

        df=pd.read_csv(f"{folder}/{LoadPath}")
        df=df[df["key"].isin(SelectedKeys)].reset_index(drop=True)

        columns = Columns if Columns is not None else _estimator_columns(df)
        for column in columns:
            if column not in df.columns:
                continue
            # PeakOptだけはフィット結果のプロットを出す（従来の挙動）。
            try:
                df = general.TempCalib(
                    df,
                    ValueKey=column,
                    ResultKey=_calibrated_name(column),
                    Title=f"{os.path.basename(folder)} / {column}",
                    plot=plot and column == "PeakOpt",
                )
            except Exception as exc:
                print(f"[警告] {column} のTempCalibに失敗しました: {exc}")

        df.to_csv(f"{folder}/{SavePath}",index=False)
        results[folder] = df

        summary = general.ResolutionSummary(
            df, columns=_comparison_columns(df, columns)
        )
        print(f"\n--- {os.path.basename(folder)}: FWHM before/after baseline calibration ---")
        print(summary.to_string(index=False))

    return results


def _comparison_columns(df, estimators):
    columns = []
    for column in estimators:
        if column in df.columns:
            columns.append(column)
        calibrated = _calibrated_name(column)
        if calibrated in df.columns:
            columns.append(calibrated)
    return columns


def CompareEstimators(csv_path, Columns=None, binNum=None, plot=True):
    """1つのCSVについて mean/std/FWHM/FWHM|mean| の比較表を返す（項目7）。

    ``Base vs Peak`` / ``Base vs PeakOpt`` / ``Base vs PeakOptTemp`` /
    ``Decay vs PeakOpt`` の散布図とヒストグラムも合わせて表示する。
    """
    df = pd.read_csv(csv_path)
    if Columns is None:
        Columns = _comparison_columns(df, _estimator_columns(df))
    Columns = [col for col in Columns if col in df.columns]
    if not Columns:
        print(f"比較できるカラムがありません: {csv_path}")
        return pd.DataFrame()

    summary = general.ResolutionSummary(df, columns=Columns, bin_num=binNum)
    print(f"\n--- {csv_path} ---")
    print(summary.to_string(index=False))

    if plot:
        DiagnosticPlots(df, Columns, binNum=binNum, title=os.path.basename(csv_path))
    return summary


def DiagnosticPlots(df, Columns=None, binNum=None, title=None):
    """相関図とヒストグラムを1枚にまとめて表示する。"""
    if Columns is None:
        Columns = _comparison_columns(df, _estimator_columns(df))
    Columns = [col for col in Columns if col in df.columns]
    if not Columns:
        return

    scatter_pairs = [("Base", col) for col in Columns]
    if "Decay" in df.columns and "PeakOpt" in df.columns:
        scatter_pairs.append(("Decay", "PeakOpt"))

    rows = 2
    cols = max(len(scatter_pairs), len(Columns))
    fig, axes = plt.subplots(rows, cols, figsize=(3.2 * cols, 6.4), squeeze=False)
    if title:
        fig.suptitle(title)

    for ax in axes.ravel():
        ax.set_visible(False)

    for i, (xkey, ykey) in enumerate(scatter_pairs):
        ax = axes[0][i]
        ax.set_visible(True)
        ax.plot(df[xkey], df[ykey], "o", markersize=1)
        ax.set_xlabel(xkey)
        ax.set_ylabel(ykey)
        ax.grid(True)

    for i, column in enumerate(Columns):
        ax = axes[1][i]
        ax.set_visible(True)
        values = df[column].to_numpy(dtype=float)
        values = values[np.isfinite(values)]
        if len(values) == 0:
            continue
        bins = binNum if binNum is not None else general.OptimalBinCount(values)
        ax.hist(values, bins=bins)
        try:
            fwhm, popt, edges = general.GaussianFWHM(values, bin_num=bins)
            x_fit = np.linspace(edges[0], edges[-1], 500)
            ax.plot(x_fit, general.gaussian(x_fit, *popt), color="red", alpha=0.6)
            ax.set_title(f"{column}\nFWHM/mean={fwhm / abs(popt[1]) * 100:.3f}%", fontsize=9)
        except Exception:
            ax.set_title(column, fontsize=9)
        ax.grid(True)

    fig.tight_layout()
    plt.show()
    plt.close(fig)


def CompareChain(config: dict, path: str, Channel, SelectedKeys=None,
                 LoadPath="output_tempcalib.csv", FilterMethod=None):
    """Getparaとの差がどの段階で出るかを順に確認する（項目8）。

    1. Average Pulse       波高・立ち上がり・raw平均との差
    2. Noise spectrum      使用中のASDの帯域とレベル
    3. Optimal Filter template  時間領域テンプレートのノルム
    4. 各イベントのOF出力  推定量ごとのmean/std
    5. Baseline補正前FWHM
    6. Baseline補正後FWHM
    """
    rate = config["Readout"]["Rate"]
    sample = config["Readout"]["Sample"]
    methods = (
        list(general.OPTIMAL_FILTER_METHODS)
        if FilterMethod is None else [FilterMethod]
    )

    print("=" * 70)
    print(f"CompareChain: CH{Channel}")
    print("=" * 70)

    # --- 1. Average pulse ---------------------------------------------------
    for method in methods:
        avg_path = OptimalFilterPath(path, Channel, method, name="average_pulse")
        raw_path = OptimalFilterPath(path, Channel, method, name="average_pulse_raw")
        if not os.path.exists(avg_path):
            continue
        avg = np.loadtxt(avg_path)
        peak_av, peak_index = general.PeakHeight(avg, config)
        print(f"\n[1] Average pulse ({method})")
        print(f"    file        : {avg_path}")
        print(f"    peak height : {peak_av:.6g} @ index {peak_index} "
              f"({peak_index / rate:.6g} s)")
        print(f"    baseline rms: {np.std(avg[:general.BaselineWindow(config, len(avg))[1]]):.6g}")
        if os.path.exists(raw_path):
            raw_avg = np.loadtxt(raw_path)
            raw_peak, _ = general.PeakHeight(raw_avg, config)
            print(f"    raw (no Bessel) peak: {raw_peak:.6g} "
                  f"(diff {100 * (peak_av - raw_peak) / raw_peak:+.3f}%)")

    # --- 2. Noise spectrum --------------------------------------------------
    for method in methods:
        noise_path = NoiseModelPath(path, Channel, method, config)
        if not os.path.exists(noise_path):
            continue
        asd = np.loadtxt(noise_path)
        fq = np.fft.rfftfreq(sample, d=1 / rate)[: len(asd)]
        print(f"\n[2] Noise ASD ({method})")
        print(f"    file   : {noise_path}")
        print(f"    length : {len(asd)} (rfft length {sample // 2 + 1})")
        print(f"    units  : ASD [pA/sqrt(Hz)]  ->  PSD = ASD**2")
        for target in (100.0, 1000.0, 10000.0):
            if target <= fq[-1]:
                idx = int(np.searchsorted(fq, target))
                idx = min(idx, len(asd) - 1)
                print(f"    ASD({fq[idx]:8.1f} Hz) = {asd[idx]:.6g}")

    # --- 3. Templates -------------------------------------------------------
    for method in methods:
        template_path = OptimalFilterPath(path, Channel, method)
        if not os.path.exists(template_path):
            continue
        template = np.loadtxt(template_path)
        print(f"\n[3] Template ({method})")
        print(f"    file  : {template_path}")
        print(f"    |h|   : {np.linalg.norm(template):.6g}, "
              f"sum={np.sum(template):.6g}, max={np.max(np.abs(template)):.6g}")

    # --- 4-6. Estimators ----------------------------------------------------
    csv_path = os.path.join(path, f"CH{Channel}_pulse", LoadPath)
    if not os.path.exists(csv_path):
        csv_path = os.path.join(path, f"CH{Channel}_pulse", "output_optimalfilter.csv")
    if not os.path.exists(csv_path):
        print(f"\n[4-6] {csv_path} が見つかりません。OptimalFilterを先に実行してください。")
        return None

    df = pd.read_csv(csv_path)
    if SelectedKeys is not None:
        df = df[df["key"].isin(SelectedKeys)].reset_index(drop=True)

    estimators = _estimator_columns(df)
    print(f"\n[4] Per-event optimal filter outputs ({len(df)} events)")
    before = general.ResolutionSummary(df, columns=estimators)
    print("\n[5] FWHM before baseline (temperature) calibration")
    print(before.to_string(index=False))

    calibrated = [_calibrated_name(col) for col in estimators]
    calibrated = [col for col in calibrated if col in df.columns]
    if calibrated:
        after = general.ResolutionSummary(df, columns=calibrated)
        print("\n[6] FWHM after baseline (temperature) calibration")
        print(after.to_string(index=False))
    else:
        after = pd.DataFrame()
        print("\n[6] 補正後カラムがありません。TempCalibを先に実行してください。")

    return {"before": before, "after": after, "csv": csv_path}

def _pulse_path(path, channel, key):
    return f"{path}/CH{channel}_pulse/rawdata/CH{channel}_{key}.dat"


def _prepared_pulse(raw, base, config, apply_bessel=True):
    """OF適用時と同じ前処理（baseline減算 → Bessel）を1か所にまとめる。

    テンプレート作成側と適用側で必ずこの関数を通すことで、両者の前処理が
    ずれないようにする。前処理が食い違うとテンプレートと波形の形が合わず、
    最適フィルタの利得が落ちてFWHMが悪化する。
    """
    pulse = np.asarray(raw, dtype=float) - float(base)
    if apply_bessel:
        cf = config["Analysis"]["CutoffFrequency"]
        if cf and cf > 0:
            pulse = general.Bessel(pulse, config["Readout"]["Rate"], cf)
    return pulse


def _resample_to(array, length):
    array = np.asarray(array, dtype=float)
    if len(array) == length:
        return array
    return np.interp(
        np.linspace(0, 1, length), np.linspace(0, 1, len(array)), array
    )


def OptimalFilter(
    config: dict,
    path: str,
    NoiseSPE,
    Channel: int,
    SelectedKeys,
    SavePath="output_optimalfilter.csv",
    FilterMethod="Current (rfft/irfft + Bessel)",
    eta_uA_per_V=None,
    Compare=True,
    plot=True,
):
    """最適フィルタを掛けて ``PeakOpt`` などの波高推定量を出力する。

    出力カラム:
      ``PeakOpt``        選択した ``FilterMethod`` の結果（本命）
      ``PeakOptLegacy``  修正前と同じ処理（raw平均パルス + テンプレートBessel）
      ``PeakOptPSD``     理論最適 ``S*/PSD`` を規格化したもの（Peakと同スケール）

    ``PeakOptLegacy`` を必ず残すのは、修正の前後を同一データ・同一
    SelectedKeys で直接比較できるようにするため（Compare=Falseで無効）。
    """
    rate = config["Readout"]["Rate"]
    sample = config["Readout"]["Sample"]

    if eta_uA_per_V is None:
        eta_uA_per_V = float(input("eta [uA/V]:"))
    else:
        eta_uA_per_V = float(eta_uA_per_V)

    df = pd.read_csv(f"{path}/CH{Channel}_pulse/output.csv")

    # --- SelectedKeys のみを残す ---
    df = df[df["key"].isin(SelectedKeys)].reset_index(drop=True)
    # 以降のキーはすべて df["key"] 由来にそろえる。SelectedKeys が float や
    # numpy 型で渡ってきても、辞書引きと df への書き戻しが食い違わない。
    keys = df["key"].to_numpy()
    base_by_key = dict(zip(keys, df["Base"].to_numpy()))

    # 平均パルスは2種類作る。
    #   average_filtered : baseline減算 → Bessel → 平均（OF適用時と同じ前処理）
    #   average_raw      : baseline減算のみ（修正前の挙動、比較用）
    average_filtered = np.zeros(sample, dtype=float)
    average_raw = np.zeros(sample, dtype=float)
    usable_keys = []

    for key in tqdm.tqdm(keys, desc="Average pulse"):
        raw = general.LoadBin(_pulse_path(path, Channel, key))
        if raw is None or len(raw) != sample:
            continue
        base = base_by_key.get(key)
        if base is None or pd.isna(base):
            print(f"[警告] key={key} に対応するBase値が見つかりません。スキップします。")
            continue
        average_raw += _prepared_pulse(raw, base, config, apply_bessel=False)
        average_filtered += _prepared_pulse(raw, base, config, apply_bessel=True)
        usable_keys.append(key)

    count = len(usable_keys)
    if count == 0:
        print("[警告] 有効なパルスがありません。")
        return
    average_raw /= count
    average_filtered /= count

    # PSD最適フィルタの出力をPeakと同じ単位で読めるようにするためのスケール。
    amplitude_scale, _peak_index = general.PeakHeight(average_filtered, config)

    templates = {
        "PeakOpt": general.OptimalFilterTemplate(
            NoiseSPE, average_filtered, config,
            method=FilterMethod, plot=plot,
            AmplitudeScale=amplitude_scale,
        )
    }
    if Compare:
        # 修正前と同じ組み合わせ: raw平均パルス + テンプレートBessel。
        templates["PeakOptLegacy"] = general.OptimalFilterTemplate(
            NoiseSPE, average_raw, config,
            method=general.CURRENT_METHOD, plot=False,
        )
        if FilterMethod != general.PSD_OPTIMAL_METHOD:
            templates["PeakOptPSD"] = general.OptimalFilterTemplate(
                NoiseSPE, average_filtered, config,
                method=general.PSD_OPTIMAL_METHOD, plot=False,
                AmplitudeScale=amplitude_scale,
            )

    templates = {name: _resample_to(filt, sample) for name, filt in templates.items()}
    filt = templates["PeakOpt"]

    # 実際に適用したテンプレートと平均パルスを保存する。
    template_path = OptimalFilterPath(path, Channel, FilterMethod)
    average_path = OptimalFilterPath(path, Channel, FilterMethod, name="average_pulse")
    np.savetxt(template_path, filt)
    np.savetxt(
        average_path,
        average_filtered,
        header="baseline-subtracted, Bessel-filtered average pulse",
    )
    np.savetxt(
        OptimalFilterPath(path, Channel, FilterMethod, name="average_pulse_raw"),
        average_raw,
        header="baseline-subtracted only (pre-fix template input)",
    )

    time = np.arange(sample) / rate
    plt.plot(time, filt)
    plt.title(f"Optimal Filter (time domain)\n{FilterMethod}")
    plt.xlabel("Time [s]")
    plt.ylabel("Amplitude")
    plt.grid()
    plt.savefig(OptimalFilterPath(path, Channel, FilterMethod, ext="png"))
    plt.cla()
    print(f"Saved optimal filter to {template_path}")
    print(f"Saved average pulse to {average_path}")

    # SelectedKeys に対応する行だけ処理
    estimates = {name: {} for name in templates}
    for key in tqdm.tqdm(usable_keys, desc="Applying filter"):
        raw = general.LoadBin(_pulse_path(path, Channel, key))
        if raw is None or len(raw) != sample:
            continue
        pulse = _prepared_pulse(raw, base_by_key[key], config, apply_bessel=True)
        for name, template in templates.items():
            estimates[name][key] = float(np.sum(pulse * template))

    for name, values in estimates.items():
        df[name] = df["key"].map(values)

    output_csv = f"{path}/CH{Channel}_pulse/{SavePath}"
    df.to_csv(output_csv, index=False)

    if Compare:
        summary = general.ResolutionSummary(df, columns=_estimator_columns(df))
        print(f"\n--- CH{Channel} pulse-height estimators (before TempCalib) ---")
        print(summary.to_string(index=False))

    return df

# ---------------------------------------------------------------------------
# Interactive views and summaries
# ---------------------------------------------------------------------------

def Scatter2D(path):
    folders = glob.glob(os.path.join(path, "CH*_pulse"))
    chs = [re.search(r'CH(.*)_pulse', os.path.basename(f)).group(1) for f in folders]

    ResultList=["Base","Peak","Rise","Decay"]

    XChannel=questionary.select("Select X Channel:",choices=chs).ask()
    XKey=questionary.select("Select X Key:",choices=ResultList).ask()

    YChannel=questionary.select("Select Y Channel:",choices=chs).ask()
    YKey=questionary.select("Select Y Key:",choices=ResultList).ask()

    dfX=pd.read_csv(f"{path}/CH{XChannel}_pulse/output.csv")
    dfY=pd.read_csv(f"{path}/CH{YChannel}_pulse/output.csv")

    dfX,dfY=general.KeyIsin(dfX,dfY)

    general.Scatter2D(dfX[XKey],dfY[YKey],xlabel=f"CH{XChannel}_{XKey}",ylabel=f"CH{YChannel}_{YKey}")

def ViewPulse(path:str,Channel:int,Key:int):
    config=general.LoadJson(f"{path}/PulseConfig.json")
    pulse=general.LoadBin(f"{path}/CH{Channel}_pulse/rawdata/CH{Channel}_{Key}.dat")
    print(f"path:{path}/CH{Channel}_pulse/rawdata/CH{Channel}_{Key}.dat")
    print(f"sample:{len(pulse)}")
    result=general.AnalyzePulse(pulse,config,Key,plot=True)
    print(result)
    
def Hist(csvpath:str,Key:str,binNum=None):
    df=pd.read_csv(csvpath)
    data=df[Key]
    fwhm,reso=general.MakeHistgram(data,bin_num=binNum)
    plt.show()
    print(f"FWHM:{fwhm}, Reso:{reso}%")


# ---------------------------------------------------------------------------
# RT analysis
# ---------------------------------------------------------------------------

def RunRT(path: str | None = None):
    if path is None:
        path = general.InputPath()
    path = os.path.abspath(path)
    if not os.path.isdir(path):
        print(f"Path not found: {path}")
        return

    with _cd(path):
        if not os.path.exists("output"):
            os.mkdir("output")
        if not os.path.exists("rawdata"):
            files = natsort.natsorted(glob.glob("*.dat"))
            os.mkdir("rawdata")
            for file_path in tqdm.tqdm(files, desc="Moving RT data", leave=False):
                shutil.move(file_path, "rawdata")
        else:
            files = natsort.natsorted(glob.glob("rawdata/*.dat"))

        I_bias = []
        V_out = []
        T = []
        for file_path in tqdm.tqdm(files, desc="Reading RT data", leave=False):
            data = _load_dat(file_path)
            name = os.path.splitext(os.path.basename(file_path))[0]
            V_out.append(np.mean(data))
            T.append(_extract_int(r"_(\d+)mK", name))
            I_bias.append(float(_extract_int(r"_(\d+)uA", name)))

        low_temp = T.count(np.min(T))

        popt, _cov = curve_fit(_linear, I_bias[:low_temp], V_out[:low_temp])
        eta_uA_per_V = 1 / popt[0]
        print(f"eta (uA/V): {eta_uA_per_V}")

        T = natsort.natsorted(set(T[low_temp:]))
        I_bias_2 = natsort.natsorted(set(I_bias[low_temp:]))

        channel_match = re.search(r"CH(\d+)_", os.path.basename(files[0])) if files else None
        channel = channel_match.group(1) if channel_match else "1"

        V_out = []
        for i in tqdm.tqdm(I_bias_2, desc="Building RT curves"):
            V = []
            for t in T:
                data = _load_dat(f"rawdata/CH{channel}_{t}mK_{int(i)}uA.dat")
                V.append(np.mean(data))
            V_out.append(V)

        V_out = np.array(V_out)
        cnt = 0
        for values in V_out:
            if cnt > 0:
                V_out_base = values - V_out[0]
                R = 3.9 * (I_bias_2[cnt] / (eta_uA_per_V * V_out_base) - 1)
                plt.title("R-T")
                plt.plot(T, R, marker="o", linewidth=1, label=f"{I_bias_2[cnt]}uA", markersize=4)
                plt.xlabel("Temperature[mK]", fontsize=16)
                plt.ylabel("Resistance[m$\\Omega$]", fontsize=16)
                plt.grid(True)
                plt.legend(loc="best", fancybox=True, shadow=True)
                np.savetxt(f"output/rt_{int(I_bias_2[cnt])}uA.txt", [T, R])
            cnt += 1

        plt.savefig("output/rt_RT.png")
        plt.show()
