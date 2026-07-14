import numpy as np
import pandas as pd
import scipy
import json
import matplotlib.pyplot as plt
import tqdm
import sys

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

def AnalyzePulse(pulse, Json: dict, key,plot=False):
    try:
        pulse = pulse.astype(float)

        pulse = Bessel(pulse, Json["Readout"]["Rate"], Json["Analysis"]["CutoffFrequency"])

        base = np.mean(pulse[0:Json["Readout"]["PreSample"]])
        pulse -= base

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

        # --- ここで有限値かチェック ---
        if not all(np.isfinite(list(result.values()))):
            return None
        if rise<0 or decay<0:
            return None
        
        if plot:
            t = np.arange(len(pulse)) / Json["Readout"]["Rate"]
            plt.figure(figsize=(10, 5))
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

def GetSelectedIndex(x, y):
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
    
    plt.plot(x, y, "bo", markersize=1)
    plt.grid()
    picked = plt.ginput(n=-1, timeout=-1)
    plt.close()

    if len(picked) < 3:
        print("3点以上選択してください")
        sys.exit()

    picked = np.array(picked)
    inside = np.zeros(len(x), dtype=bool)
    
    for i, (sx, sy) in enumerate(zip(x, y)):
        inside[i] = inpolygon(sx, sy, picked[:, 0], picked[:, 1])

    selected_index = np.where(inside)[0]
    return selected_index

def SelectIDFrom2DF(dfX,dfY,key:str):
    selected_index = GetSelectedIndex(dfX[key], dfY[key])
    selected_ids = dfX.iloc[selected_index]["key"].values
    return selected_ids

def SelectIDFrom1DF(df,keyX:str,keyY:str):
    selected_index = GetSelectedIndex(df[keyX], df[keyY])
    selected_ids = df.iloc[selected_index]["key"].values
    return selected_ids

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

def OptimalFilterTemplate(NoiseSPE, AveragePulse, config):
    rate = config["Readout"]["Rate"]
    sample = config["Readout"]["Sample"]
    cf = config["Analysis"]["CutoffFrequency"]

    # --- FFT関連の基本設定 ---
    fq = np.fft.rfftfreq(sample, d=1/rate)  # 片側スペクトルの周波数軸
    F = np.fft.rfft(AveragePulse)            # 片側FFT（実信号用）

    # --- NoiseSPE の長さを調整 ---
    if len(NoiseSPE) != len(F):
        N_target = len(F)
        N_source = len(NoiseSPE)
        x_source = np.linspace(0, 1, N_source)
        x_target = np.linspace(0, 1, N_target)
        NoiseSPE = np.interp(x_target, x_source, NoiseSPE)
        plt.loglog(np.abs(NoiseSPE))
        plt.title("Interpolated Noise Spectrum")
        plt.show()

    # --- ローパス処理 ---
    F_filtered = np.copy(F)
    F_filtered[fq > cf] = 0  # ナイキスト考慮済み

    # --- 最適フィルタ作成 ---
    filt_fft = F_filtered / NoiseSPE
    filt_time = np.fft.irfft(filt_fft, n=sample).real  # 逆FFT（片側）
    filt_time=Bessel(filt_time,rate,cf)

    # --- 時間軸を生成 ---
    time = np.arange(sample) / rate

    # --- 可視化 ---
    plt.plot(time, filt_time)
    plt.title("Optimal Filter (time domain)")
    plt.xlabel("Time [s]")
    plt.ylabel("Amplitude")
    plt.show()

    return filt_time

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

def GN(AMpModel):
    random_phases = np.exp(1j * np.random.uniform(0, 2*np.pi, len(AMpModel)))
    spec = AMpModel * random_phases
    noise= np.fft.irfft(spec)
    return noise