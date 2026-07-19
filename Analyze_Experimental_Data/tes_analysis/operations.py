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
    return filter_method == "Legacy (fft/ifft)"


def _legacy_output_name(config):
    """Return the filename used for the legacy numerical result."""
    return "modelnoise_old.txt"


def NoiseModelPath(path, channel, filter_method="Current (rfft/irfft + Bessel)", config=None):
    """Return the modelnoise path corresponding to the selected FFT policy."""
    noise_folder = os.path.join(path, f"CH{channel}_noise")
    if _is_legacy_filter_method(filter_method):
        return os.path.join(noise_folder, _legacy_output_name(config))
    return os.path.join(noise_folder, "modelnoise.txt")


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
):
    
    # 設定値の取得
    sample = config["Readout"]["Sample"] # データ点数 (N)
    rate = config["Readout"]["Rate"]     # サンプリングレート (Fs)
    cutoff = config["Analysis"]["CutoffFrequency"] # ローパスフィルタのカットオフ周波数

    eta=float(input("eta:"))

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
            np.savetxt(noise_model_path, amp_dens * eta * 1e6)
            plt.plot(
                general.GetFreq(rate, sample)[: sample // 2 + 1],
                amp_dens * eta * 1e6,
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
        
        amplitude_model = np.zeros(sample)
        count = 0
        original_noise_list = []
        
        # メディアンフィルタのカーネルサイズ (窓サイズ) を定義
        # 奇数に設定し、ノイズの幅に応じて調整します。3, 5, 7 などが一般的です。
        # ここでは例としてカーネルサイズ3を使用します。
        median_kernel_size = 3

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
            noise = general.Bessel(noise, rate, cutoff)
            diff=np.max(noise)-np.min(noise)
            if diff>noise_threshold:
                filtered+=1
                continue
            original_noise_list.append(noise)

            noise_fft = np.fft.fft(noise)
            noise_amp = np.abs(noise_fft)
            amplitude_model += noise_amp
            count += 1

        print(f"count;{count}")

        print(f"Filtered {filtered} and skipped {sample_unmatch} from {len(noise_pathes)}.")

        if count == 0:
            print(f"No usable noise data found in {folder}; skipping.")
            continue

        amplitude_model/=count

        np.savetxt(f"{folder}/noise_fft_Amplitude.txt",amplitude_model)

        df=rate/sample
        power=amplitude_model**2 / df
        amp_dens=np.sqrt(power)
        amp_dens = amp_dens[: int(sample / 2) + 1] * eta * 1e+6

        noise_model_path = NoiseModelPath(
            path,
            os.path.basename(folder).removeprefix("CH").removesuffix("_noise"),
            FilterMethod,
            config,
        )
        os.makedirs(os.path.dirname(noise_model_path), exist_ok=True)
        np.savetxt(noise_model_path, amp_dens)

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

def TempCalib(path:str,SelectedKeys,SavePath="output_tempcalib.csv",LoadPath="output_optimalfilter.csv"):
    for folder in tqdm.tqdm(glob.glob(f"{path}/CH*_pulse")):
        if os.path.exists(f"{folder}/{LoadPath}"):
            df=pd.read_csv(f"{folder}/{LoadPath}")
            df=df[df["key"].isin(SelectedKeys)].reset_index(drop=True)
            data_calib=general.TempCalib(df)
            data_calib.to_csv(f"{folder}/{SavePath}",index=False)
        else:
            print(f"output.csv not found in {folder}/{LoadPath}, skipping TempCalib.")

def OptimalFilter(
    config: dict,
    path: str,
    NoiseSPE,
    Channel: int,
    SelectedKeys,
    SavePath="output_optimalfilter.csv",
    FilterMethod="Current (rfft/irfft + Bessel)",
):
    rate = config["Readout"]["Rate"]
    cf = config["Analysis"]["CutoffFrequency"]

    eta = float(input("eta:"))

    df = pd.read_csv(f"{path}/CH{Channel}_pulse/output.csv")

    # --- SelectedKeys のみを残す ---
    df = df[df["key"].isin(SelectedKeys)].reset_index(drop=True)

    AveragePulse = np.zeros(config["Readout"]["Sample"], dtype=float)
    count = 0

    for key in tqdm.tqdm(SelectedKeys):
        pulse = general.LoadBin(f"{path}/CH{Channel}_pulse/rawdata/CH{Channel}_{key}.dat")
        if len(pulse) != config["Readout"]["Sample"]:
            continue
        AveragePulse += pulse
        count += 1

    if count > 0:
        AveragePulse /= count
    else:
        print("[警告] 有効なパルスがありません。")
        return

    filt = general.OptimalFilterTemplate(
        NoiseSPE, AveragePulse, config, method=FilterMethod
    )

    pulse_length = len(AveragePulse)
    # 各配列の「もとのインデックス軸」
    x_avg = np.linspace(0, 1, pulse_length)

    # filt の補間処理
    if len(filt) != pulse_length:
        x_filt = np.linspace(0, 1, len(filt))
        filt = np.interp(x_avg, x_filt, filt)

    # SelectedKeys に対応する行だけ処理
    for key in tqdm.tqdm(SelectedKeys):
        pulse = general.LoadBin(f"{path}/CH{Channel}_pulse/rawdata/CH{Channel}_{key}.dat").copy()
        if len(pulse) != config["Readout"]["Sample"]:
            continue
        base_values = df.loc[df["key"] == key, "Base"].values

        if len(base_values) > 0 and not pd.isna(base_values[0]):
            pulse -= base_values[0]
        else:
            print(f"[警告] key={key} に対応するBase値が見つかりません。スキップします。")
            continue

        pulse = general.Bessel(pulse, rate, cf)
        df.loc[df["key"] == key, "PeakOpt"] = np.sum(pulse * filt)

    df.to_csv(f"{path}/CH{Channel}_pulse/{SavePath}", index=False)

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
        eta = 1 / popt[0]
        print(f"eta (uA/V): {eta}")

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
                R = 3.9 * (I_bias_2[cnt] / (eta * V_out_base) - 1)
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
