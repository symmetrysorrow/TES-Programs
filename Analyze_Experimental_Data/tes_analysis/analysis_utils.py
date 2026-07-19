"""Shared helpers for TES data analysis.

The legacy scripts import this module directly, so its public function names
are kept stable.  Functions are grouped by responsibility below: file I/O,
signal processing, pulse analysis, interactive selection/plotting, and noise
synthesis.
"""

import json

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy
import tqdm


# ---------------------------------------------------------------------------
# File I/O and common axes
# ---------------------------------------------------------------------------

def LoadTxt(file_path:str):
    try:
        data = np.loadtxt(file_path, comments="#")
        return data
    except Exception as e:
        print(f"Error loading file {file_path}: {e}")
        return None
    
def LoadBin(file_path:str):
    try:
        with open(file_path, "rb") as fb:
            fb.seek(4)
            data = np.frombuffer(fb.read(), dtype="float64")
        return data
    except Exception:
        try:
            data=LoadTxt(file_path)
            return data
        except Exception as e:
            print(f"Error loading binary file {file_path}: {e}")
            return None
    
def LoadJson(file_path:str):
    try:
        with open(file_path, "r") as f:
            data = json.load(f)
        return data
    except Exception as e:
        print(f"Error loading JSON file {file_path}: {e}")
        return None
    
def Bessel(data,rate:float,fs:float):
    ws=GetWs(rate,fs)
    b,a=scipy.signal.bessel(2,ws,"low")
    filtered_data=scipy.signal.filtfilt(b,a,data)
    return filtered_data

def gaussian(x, amp, mean, stddev):
    return amp * np.exp(-((x - mean) ** 2) / (2 * stddev ** 2))

def OptimalBinCount(data):
    q1, q3 = np.percentile(data, [25, 75])
    iqr = q3 - q1  # 四分位範囲
    bin_width = 2 * iqr / (len(data) ** (1/3))  # ビン幅
    bin_count = int(np.ceil((np.max(data) - np.min(data)) / bin_width))  # ビン数
    return max(bin_count, 1)  # ビン数が1未満にならないようにする

def MakeHistgram(data,bin_num=None,label=None,HistColor=None):
    if bin_num is None:
        bin_num = OptimalBinCount(data)
    hist, bin_edges = np.histogram(data, bins=bin_num, density=False)
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2  # ビンの中心を計算
    initial_guess = [np.max(hist), np.mean(data), np.std(data)]
    if HistColor is not None:
        plt.hist(data, bins=bin_num, density=False, label=label,color=HistColor)
    else:
        plt.hist(data, bins=bin_num, density=False, label=label)
    # ガウスフィッティング
    popt, pcov = scipy.optimize.curve_fit(gaussian, bin_centers, hist, p0=initial_guess, maxfev=10000000)
    amp_fit, mean_fit, stddev_fit = popt
    fwhm = 2 * stddev_fit * np.sqrt(2 * np.log(2))

    #ヒストグラム
    x_fit = np.linspace(bin_edges[0], bin_edges[-1], 1000)  # フィッティング用のx
    plt.plot(x_fit, gaussian(x_fit, *popt),color="red",alpha=0.5)  # フィッティング曲線

    return fwhm,fwhm/mean_fit

def InputPath():
    path=input("Input file path:")
    return path

def GetWs(rate:float,fs:float):
    ws=fs/rate*2
    return ws

def GetFreq(rate:float,samples:int, Nyquist=False):
    fq = np.arange(0, rate, rate / samples)
    return fq

def GetTime(rate:float,samples:int):
    time = np.arange(0, samples) / rate
    return time


# ---------------------------------------------------------------------------
# Pulse analysis and calibration
# ---------------------------------------------------------------------------

def AnalyzePulse(pulse, Json: dict, key,plot=False):
    try:
        pulse = pulse.astype(float)      

        base = np.mean(pulse[0:Json["Readout"]["PreSample"]])
        pulse -= base

        rawpulse=pulse.copy()

        pulse = Bessel(pulse, Json["Readout"]["Rate"], Json["Analysis"]["CutoffFrequency"])

        peak_index = np.argmax(pulse)
        peak_av = np.mean(
            pulse[
                peak_index - Json["Analysis"]["PeakAveragePreSample"] : 
                peak_index + Json["Analysis"]["PeakAveragePostSample"]
            ]
        )

        rise_high = rise_low = 0
        for i in reversed(range(0, peak_index)):
            if pulse[i] <= peak_av * Json["Analysis"]["RiseHighRatio"]:
                rise_high = i
                break
        for j in reversed(range(0, rise_high)):
            if pulse[j] <= peak_av * Json["Analysis"]["RiseLowRatio"]:
                rise_low = j
                break
        rise = (rise_high - rise_low) / Json["Readout"]["Rate"]

        decay_high = decay_low = 0
        for i in range(peak_index, len(pulse)):
            if pulse[i] <= peak_av * Json["Analysis"]["DecayHighRatio"]:
                decay_high = i
                break
        for j in range(decay_high, len(pulse)):
            if pulse[j] <= peak_av * Json["Analysis"]["DecayLowRatio"]:
                decay_low = j
                break
        decay = (decay_low - decay_high) / Json["Readout"]["Rate"]

        result = {
            "key": int(key),
            "Base": float(base),
            "Peak": float(peak_av),
            "Rise": float(rise),
            "Decay": float(decay),
        }
 
        if plot:
            t = np.arange(len(pulse)) / Json["Readout"]["Rate"]
            plt.figure(figsize=(10, 5))
            plt.plot(t, rawpulse, label="Raw Pulse", color="lightgray")
            plt.plot(t, pulse, label="Pulse", color="gray")

            # --- 範囲を安全に切り詰める ---
            rise_low = max(0, min(len(pulse) - 1, rise_low))
            rise_high = max(0, min(len(pulse) - 1, rise_high))
            decay_low = max(0, min(len(pulse) - 1, decay_low))
            decay_high = max(0, min(len(pulse) - 1, decay_high))
            peak_pre = max(0, peak_index - Json["Analysis"]["PeakAveragePreSample"])
            peak_post = min(len(pulse) - 1, peak_index + Json["Analysis"]["PeakAveragePostSample"])

            # --- 範囲ハイライト ---
            plt.axvspan(t[rise_low], t[rise_high], color="lime", alpha=0.3, label="Rise")
            plt.axvspan(t[peak_pre], t[peak_post], color="orange", alpha=0.3, label="Peak")
            plt.axvspan(t[decay_high], t[decay_low], color="deepskyblue", alpha=0.3, label="Decay")

            # --- 代表点マーク ---
            plt.scatter([t[peak_index]], [pulse[peak_index]], color="red", marker="^", zorder=3,label="Peak")
            peak_center_time = (t[peak_pre] + t[peak_post]) / 2
            plt.scatter([peak_center_time], [peak_av], color="darkorange", marker="D", zorder=4, label="PeakAverage")

            plt.title(f"Pulse {key} - Rise/Decay Analysis")
            plt.xlabel("Time [s]")
            plt.ylabel("Amplitude")
            plt.legend()
            plt.grid(True)
            plt.tight_layout()
            plt.show()

        # --- ここで有限値かチェック ---
        if not all(np.isfinite(list(result.values()))):
            return None
        if rise<0 or decay<0:
            return None

        return result

    except Exception as e:
        print(f"Error in AnalyzePulse: {e}")
        return None

def TempCalib(data):
    def func(X, *params):
        Y = np.zeros_like(X)
        for i, param in enumerate(params):
            Y = Y + np.array(param * X ** i)
        return Y
    
    def Calibration(x,params):
        array = np.zeros(len(params))
        for i,param in enumerate(params):
            term = param * x ** i
            array[i] = term
            sum = np.sum(array)
        return sum
    
    bases=data["Base"]
    heights_opt=data["PeakOpt"]

    p0=[0.01,0.01,0.01,0.01,0.01,0.01]

    popt,_pcov = scipy.optimize.curve_fit(func,bases,heights_opt,p0)
    x_fit = np.linspace(np.min(bases),np.max(bases),100000)
    fitted = func(x_fit,*tuple(popt))

    plt.plot(bases,heights_opt,'o',color='blue',markersize=3,label='a')
    plt.plot(x_fit,fitted,color='red',linewidth=1.0,linestyle='-')
    plt.xlabel('baseline [V]',fontsize = 16)
    plt.ylabel('pulseheight [V]',fontsize = 16)
    plt.grid()
    plt.show()
    plt.cla()

    st=np.mean(heights_opt)

    for index,row in tqdm.tqdm(data.iterrows()):
        data.at[index,"PeakOptTemp"] = row['Peak']/Calibration(row['Base'],popt)*st

    plt.plot(bases,data["PeakOptTemp"],'o',color='tab:blue',markersize=0.7,label='a')
    plt.xlabel('baseline [V]',fontsize = 16)
    plt.ylabel('pulseheight [V]',fontsize = 16)
    plt.grid()
    plt.show()
    plt.cla()

    return data


# ---------------------------------------------------------------------------
# Interactive selection and plotting
# ---------------------------------------------------------------------------

def GetSelectedKey(xDf, yDf, xkey, ykey, key_col="key"):
    def inpolygon(sx, sy, poly_x, poly_y):
        inside = False
        n = len(poly_x)
        j = n - 1
        for i in range(n):
            xi, yi = poly_x[i], poly_y[i]
            xj, yj = poly_x[j], poly_y[j]
            if ((yi > sy) != (yj > sy)) and \
               (sx < (xj - xi) * (sy - yi) / (yj - yi) + xi):
                inside = not inside
            j = i
        return inside

    # --- 共通キーのみに絞る ---
    common_keys = np.intersect1d(xDf[key_col].values, yDf[key_col].values)
    xDf = xDf[xDf[key_col].isin(common_keys)].reset_index(drop=True)
    yDf = yDf[yDf[key_col].isin(common_keys)].reset_index(drop=True)

    if len(xDf) == 0:
        print("共通するキーが存在しません。")
        return np.array([])

    # --- x, y, key 抽出 ---
    x = xDf[xkey].values
    y = yDf[ykey].values
    keys = xDf[key_col].values

    if len(x) != len(y):
        raise ValueError(f"x({len(x)}) と y({len(y)}) の長さが一致していません")

    # --- プロットと選択 ---
    fig, ax = plt.subplots()
    ax.plot(x, y, "bo", markersize=3)
    ax.set_xlabel(xkey)
    ax.set_ylabel(ykey)
    ax.grid()

    picked = []

    def onclick(event):
        toolbar = plt.get_current_fig_manager().toolbar
        if toolbar.mode != '':  # ズーム・パン中は無視
            return
        if event.inaxes != ax:
            return
        picked.append((event.xdata, event.ydata))
        ax.plot(event.xdata, event.ydata, "r+", markersize=8)
        fig.canvas.draw()

    def onclose(event):
        print("ウィンドウが閉じられました。選択を終了します。")

    fig.canvas.mpl_connect('button_press_event', onclick)
    fig.canvas.mpl_connect('close_event', onclose)

    print("クリックで点を選択（ズーム・パン中は無視）")
    print("ウィンドウを閉じると選択終了")

    plt.show()  # GUIループ。閉じると自動的に続行される

    # --- 選択領域判定 ---
    if len(picked) < 3:
        print("3点以上選択してください")
        return np.array([])

    picked = np.array(picked)
    inside = np.zeros(len(x), dtype=bool)

    for i, (sx, sy) in enumerate(zip(x, y)):
        inside[i] = inpolygon(sx, sy, picked[:, 0], picked[:, 1])

    selected_keys = xDf.loc[inside, key_col].values

    #print(f"Selected keys ({len(selected_keys)}件): {selected_keys}")
    return selected_keys


def SelectIDFrom2DF(dfX,dfY,key:str):
    selected_key = GetSelectedKey(dfX, dfY, key, key)
    #selected_ids = dfX.iloc[selected_index]["key"].values
    #return selected_ids
    return selected_key

def SelectIDFrom1DF(df,keyX:str,keyY:str):
    # GetSelectedKey already returns values from the ``key`` column.  Treating
    # those key values as positional row indices makes ``iloc`` fail whenever
    # a selected key is outside the DataFrame's positional range.
    return GetSelectedKey(df, df, keyX, keyY)

def Scatter2D(x,y,xlabel=None,ylabel=None,title=None):
    plt.plot(x, y, "bo", markersize=1)
    if xlabel is not None:
        plt.xlabel(xlabel)
    if ylabel is not None:
        plt.ylabel(ylabel)
    if title is not None:
        plt.title(title)
    plt.grid()
    plt.show()
    plt.cla()

def KeyIsin(df1,df2):
    common_keys = set(df1["key"]) & set(df2["key"])

    df1=df1[df1["key"].isin(common_keys)].reset_index(drop=True)
    df2=df2[df2["key"].isin(common_keys)].reset_index(drop=True)
    return df1,df2


# ---------------------------------------------------------------------------
# Optimal filtering
# ---------------------------------------------------------------------------

def _resize_spectrum(spectrum, target_length):
    """Resize a real-valued spectrum while keeping the existing input format."""
    spectrum = np.asarray(spectrum)
    if len(spectrum) == target_length:
        return spectrum
    x_source = np.linspace(0, 1, len(spectrum))
    x_target = np.linspace(0, 1, target_length)
    return np.interp(x_target, x_source, spectrum)


def OptimalFilterTemplate(
    NoiseSPE,
    AveragePulse,
    config,
    method="Current (rfft/irfft + Bessel)",
):
    rate = config["Readout"]["Rate"]
    sample = config["Readout"]["Sample"]
    cf = config["Analysis"]["CutoffFrequency"]

    if method == "Legacy (fft/ifft)":
        # 旧版と同じ全周波数FFT・旧周波数軸・ifftの組み合わせ。
        # 現行のmodelnoise.txtは片側長の場合があるため、除算できる全長へ
        # 補間する（入力ファイル形式自体は変更しない）。
        fq = np.arange(0, rate, rate / sample)
        F = scipy.fftpack.fft(AveragePulse)
        NoiseSPE = _resize_spectrum(NoiseSPE, sample)
        F_filtered = np.copy(F)
        F_filtered[fq > cf] = 0
        filt_time = scipy.fftpack.ifft(F_filtered / NoiseSPE).real
    else:
        if method != "Current (rfft/irfft + Bessel)":
            raise ValueError(f"Unknown optimal filter method: {method}")

        # 現行方式：実信号用の片側FFTとirfftを使用。
        fq = np.fft.rfftfreq(sample, d=1 / rate)
        F = np.fft.rfft(AveragePulse)
        NoiseSPE = _resize_spectrum(NoiseSPE, len(F))
        F_filtered = np.copy(F)
        F_filtered[fq > cf] = 0
        filt_time = np.fft.irfft(F_filtered / NoiseSPE, n=sample).real
        filt_time = Bessel(filt_time, rate, cf)

    # --- 時間軸を生成 ---
    time = np.arange(sample) / rate

    # --- 可視化 ---
    plt.plot(time, filt_time)
    plt.title("Optimal Filter (time domain)")
    plt.xlabel("Time [s]")
    plt.ylabel("Amplitude")
    plt.show()

    return filt_time


# ---------------------------------------------------------------------------
# Noise synthesis
# ---------------------------------------------------------------------------

def RandomNoise(noise_fft,rate):
    random_phase = np.random.uniform(0, 2*np.pi, int(rate/2)+1)
    # 片側スペクトル（DCとNyquistは実数）
    X_half = noise_fft[:int(rate/2)+1] * np.exp(1j * random_phase)
    X_half[0] = noise_fft[0]  # DC
    if rate % 2 == 0:
        X_half[-1] = noise_fft[int(rate/2)]  # Nyquist

    # 両側スペクトルを構築（共役対称）
    X_full = np.zeros(rate, dtype=complex)
    X_full[:int(rate/2)+1] = X_half
    X_full[int(rate/2)+1:] = np.conj(X_half[-2:0:-1])

    # ifftによる時間波形の再構成
    noise_reconstructed = np.fft.ifft(X_full).real

    return noise_reconstructed

def RandomNoiseN(noise_fft):
    rate=len(noise_fft)
    random_phase = np.random.uniform(0, 2*np.pi, int(rate/2)+1)
    # 片側スペクトル（DCとNyquistは実数）
    X_half = noise_fft[:int(rate/2)+1] * np.exp(1j * random_phase)
    X_half[0] = noise_fft[0]  # DC
    if rate % 2 == 0:
        X_half[-1] = noise_fft[int(rate/2)]  # Nyquist

    # 両側スペクトルを構築（共役対称）
    X_full = np.zeros(rate, dtype=complex)
    X_full[:int(rate/2)+1] = X_half
    X_full[int(rate/2)+1:] = np.conj(X_half[-2:0:-1])

    # ifftによる時間波形の再構成
    noise_reconstructed = np.fft.ifft(X_full).real

    return noise_reconstructed

import numpy as np

def GenerateNoiseFromModel(amp_dens, sample: int, rate: float, eta: float = 1.0) -> np.ndarray:
    # --- 1. モデルノイズ読み込み ---
    n_half = len(amp_dens)
    
    # --- 2. 周波数分解能を求める ---
    df = rate / sample
    
    # --- 3. パワースペクトル密度から振幅スペクトルに戻す ---
    # amp_dens = sqrt(power) = amplitude_model / sqrt(df)
    power = (amp_dens / (eta * 1e6))**2  # もとのpower [V^2/Hz]に戻す
    amplitude_model = np.sqrt(power * df)

    # --- 4. ランダム位相の生成 ---
    random_phases = np.exp(1j * np.random.uniform(0, 2*np.pi, n_half))

    # --- 5. 片側スペクトルから両側へ ---
    spec_half = amplitude_model * random_phases
    spec_full = np.zeros(sample, dtype=complex)
    spec_full[:n_half] = spec_half
    # 実信号にするため共役対称化
    spec_full[n_half:] = np.conj(spec_half[-2:0:-1])

    # --- 6. 時間領域に変換 ---
    noise = np.fft.ifft(spec_full).real

    # --- 7. 正規化（任意） ---
    noise -= np.mean(noise)
    noise /= np.std(noise)

    return noise

def GN(AMpModel):
    random_phases = np.exp(1j * np.random.uniform(0, 2*np.pi, len(AMpModel)))
    spec = AMpModel * random_phases
    noise = np.fft.irfft(spec)
    return noise
