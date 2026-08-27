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
    
# 点数として使う設定値（np.zerosやスライスに渡すのでintでなければならない）。
# PulseConfig.jsonに 50000.0 のような実数で書かれていても動くようにする。
_SAMPLE_COUNT_FIELDS = {
    "Readout": ("Sample", "PreSample"),
    "Analysis": (
        "BaseLinePreSample",
        "BaseLinePostSample",
        "PeakAveragePreSample",
        "PeakAveragePostSample",
        "PeakSearchSample",
        "BaseStart",
        "BaseWidth",
    ),
}


def NormalizeConfig(config):
    """点数系の設定値をintへ揃える。"""
    if not isinstance(config, dict):
        return config
    for section, fields in _SAMPLE_COUNT_FIELDS.items():
        block = config.get(section)
        if not isinstance(block, dict):
            continue
        for field in fields:
            value = block.get(field)
            if isinstance(value, float) and value.is_integer():
                block[field] = int(value)
    return config


def LoadJson(file_path:str):
    try:
        with open(file_path, "r") as f:
            data = json.load(f)
        return NormalizeConfig(data)
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
    # 分布が縮退している（IQR=0 や全点同値）と0除算/NaNになるので既定値へ逃がす。
    span = np.max(data) - np.min(data)
    if not np.isfinite(bin_width) or bin_width <= 0 or not np.isfinite(span):
        return max(int(np.sqrt(len(data))), 1)
    bin_count = int(np.ceil(span / bin_width))  # ビン数
    return max(bin_count, 1)  # ビン数が1未満にならないようにする

def GaussianFWHM(data, bin_num=None):
    """ヒストグラムをガウスでフィットして ``(fwhm, popt, bin_edges)`` を返す。

    プロットを伴わない純関数なので、分解能の一括比較から呼べる。
    """
    data = np.asarray(data, dtype=float)
    data = data[np.isfinite(data)]
    if len(data) < 3:
        raise ValueError("not enough finite samples for a Gaussian fit")
    if bin_num is None:
        bin_num = OptimalBinCount(data)
    hist, bin_edges = np.histogram(data, bins=bin_num, density=False)
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    initial_guess = [np.max(hist), np.mean(data), np.std(data)]
    popt, _pcov = scipy.optimize.curve_fit(
        gaussian, bin_centers, hist, p0=initial_guess, maxfev=10000000
    )
    fwhm = 2 * abs(popt[2]) * np.sqrt(2 * np.log(2))
    return fwhm, popt, bin_edges


def MakeHistgram(data,bin_num=None,label=None,HistColor=None):
    if bin_num is None:
        bin_num = OptimalBinCount(data)
    if HistColor is not None:
        plt.hist(data, bins=bin_num, density=False, label=label,color=HistColor)
    else:
        plt.hist(data, bins=bin_num, density=False, label=label)
    # ガウスフィッティング
    fwhm, popt, bin_edges = GaussianFWHM(data, bin_num=bin_num)
    mean_fit = popt[1]

    #ヒストグラム
    x_fit = np.linspace(bin_edges[0], bin_edges[-1], 1000)  # フィッティング用のx
    plt.plot(x_fit, gaussian(x_fit, *popt),color="red",alpha=0.5)  # フィッティング曲線

    return fwhm,fwhm/mean_fit


def ResolutionSummary(data, columns=None, bin_num=None):
    """列ごとに mean / std / FWHM / FWHM|mean| をまとめた DataFrame を返す。

    ``data`` は DataFrame でも ``{name: array}`` でもよい。ガウスフィットに
    失敗した列は FWHM を NaN にして残す（比較表から列が消えない方が良い）。
    """
    if isinstance(data, pd.DataFrame):
        frame = data
        if columns is None:
            columns = [
                col for col in frame.columns
                if col != "key" and pd.api.types.is_numeric_dtype(frame[col])
            ]
        series = {col: frame[col].to_numpy(dtype=float) for col in columns if col in frame}
    else:
        series = {name: np.asarray(values, dtype=float) for name, values in data.items()}
        if columns is not None:
            series = {name: series[name] for name in columns if name in series}

    rows = []
    for name, values in series.items():
        values = values[np.isfinite(values)]
        if len(values) == 0:
            continue
        mean = float(np.mean(values))
        try:
            fwhm, popt, _edges = GaussianFWHM(values, bin_num=bin_num)
            fit_mean = float(popt[1])
        except Exception:
            fwhm = np.nan
            fit_mean = np.nan
        ratio = fwhm / abs(fit_mean) if fit_mean not in (0.0,) and np.isfinite(fit_mean) else np.nan
        rows.append(
            {
                "column": name,
                "N": int(len(values)),
                "mean": mean,
                "std": float(np.std(values, ddof=1)) if len(values) > 1 else np.nan,
                "fit_mean": fit_mean,
                "FWHM": fwhm,
                "FWHM/mean": ratio,
                "Reso[%]": ratio * 100.0 if np.isfinite(ratio) else np.nan,
            }
        )
    return pd.DataFrame(rows)

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

def BaselineWindow(config: dict, length=None):
    """ベースライン平均に使う ``[start, stop)`` を返す。

    既定（``BaseStart``/``BaseWidth`` 未設定）は従来どおり pre-trigger 全体
    ``pulse[0:PreSample]``。Getpara は ``data[presamples-base_x :
    presamples-base_x+base_w]`` を使うので、そちらへ合わせたいときは
    ``Analysis.BaseStart`` (=base_x) と ``Analysis.BaseWidth`` (=base_w) を
    設定する。例: presamples=1000 なら ``BaseStart=1000, BaseWidth=500``。

    ベースライン区間はゲイン補正 (TempCalib) の説明変数そのものなので、
    区間が長すぎる／パルスの立ち上がりを含むと Base のばらつきが増え、
    補正後の分解能が悪化する。
    """
    readout = config.get("Readout", {}) if isinstance(config, dict) else {}
    analysis = config.get("Analysis", {}) if isinstance(config, dict) else {}

    presample = int(readout.get("PreSample", 0) or 0)
    base_start = analysis.get("BaseStart")
    base_width = analysis.get("BaseWidth")

    if base_start is None or base_width is None:
        start, stop = 0, presample
    else:
        start = presample - int(base_start)
        stop = start + int(base_width)

    limit = presample if length is None else int(length)
    start = max(0, min(start, limit))
    stop = max(start + 1, min(stop, limit))
    return start, stop


def ComputeBaseline(pulse, config: dict):
    """``BaselineWindow`` で決まる区間の平均を返す。

    テンプレート作成側と適用側で必ず同じ定義を使うため、baseline は常に
    この関数を経由する。
    """
    pulse = np.asarray(pulse, dtype=float)
    start, stop = BaselineWindow(config, len(pulse))
    return float(np.mean(pulse[start:stop]))


def PeakSearchWindow(config: dict, length):
    """波高探索に使う ``[start, stop)`` を返す。

    ``Analysis.PeakSearchSample`` が正なら Getpara と同じく trigger 以降
    ``[PreSample, PreSample+PeakSearchSample)`` に限定する。全区間 argmax は
    pre-trigger のノイズスパイクや後続イベントを拾って波高を壊すことがある。
    """
    readout = config.get("Readout", {}) if isinstance(config, dict) else {}
    analysis = config.get("Analysis", {}) if isinstance(config, dict) else {}

    length = int(length)
    width = analysis.get("PeakSearchSample")
    if not width or int(width) <= 0:
        return 0, length

    start = int(readout.get("PreSample", 0) or 0)
    start = max(0, min(start, length - 1))
    stop = min(length, start + int(width))
    if stop <= start:
        return 0, length
    return start, stop


def PeakHeight(pulse, config: dict):
    """フィルタ済み・baseline減算済み波形から (peak_av, peak_index) を返す。"""
    pulse = np.asarray(pulse, dtype=float)
    analysis = config.get("Analysis", {})
    start, stop = PeakSearchWindow(config, len(pulse))
    peak_index = int(np.argmax(pulse[start:stop])) + start

    pre = int(analysis.get("PeakAveragePreSample", 0) or 0)
    post = int(analysis.get("PeakAveragePostSample", 0) or 0)
    lo = max(0, peak_index - pre)
    hi = min(len(pulse), peak_index + post)
    if hi <= lo:
        hi = lo + 1
    return float(np.mean(pulse[lo:hi])), peak_index


def AnalyzePulse(pulse, Json: dict, key,plot=False):
    try:
        pulse = pulse.astype(float)

        base = ComputeBaseline(pulse, Json)
        pulse -= base

        rawpulse=pulse.copy()

        pulse = Bessel(pulse, Json["Readout"]["Rate"], Json["Analysis"]["CutoffFrequency"])

        peak_av, peak_index = PeakHeight(pulse, Json)

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

def _polynomial(X, *params):
    """``sum(params[i] * X**i)``。curve_fit と補正の両方で使う唯一の定義。"""
    X = np.asarray(X, dtype=float)
    Y = np.zeros_like(X)
    for i, param in enumerate(params):
        Y = Y + param * X ** i
    return Y


def TempCalib(data, ValueKey="PeakOpt", ResultKey="PeakOptTemp", Title=None, plot=True):
    """ベースライン依存のゲインを1次直線でフィットして補正する。

    高次多項式だと点が疎な領域でフィット曲線が波打ち、補正後に不自然な
    構造を作るため、直線に固定している。

    フィットに使う量 (``Base`` vs ``ValueKey``) と、補正で割る量は必ず同じ
    ``ValueKey`` でなければならない。Getpara の ``temp_calib.py`` は
    ``height_opt`` でフィットして ``height_opt`` を補正しており、ここで
    ``Peak`` のような別カラムを割ると補正が意味を成さず分解能が壊れる。
    """
    # ValueKeyとResultKeyが同じ（繰り返し補正）でも壊れないよう値を控えておく。
    bases = data["Base"].to_numpy(dtype=float, copy=True)
    heights_opt = data[ValueKey].to_numpy(dtype=float, copy=True)

    finite = np.isfinite(bases) & np.isfinite(heights_opt)
    if np.count_nonzero(finite) < 2:
        raise ValueError(f"TempCalib: {ValueKey} に有効な値がありません")

    # 1次直線（切片と傾き）。パラメータを増やすとフィットがぐにゃぐにゃになる。
    p0 = [0.01, 0.01]

    popt, _pcov = scipy.optimize.curve_fit(
        _polynomial, bases[finite], heights_opt[finite], p0
    )

    st = np.mean(heights_opt[finite])
    # Getpara と同じ height/f(base)*mean(height)。ループではなくベクトル化。
    gain = _polynomial(bases, *popt)
    data[ResultKey] = heights_opt / gain * st

    if plot:
        x_fit = np.linspace(np.min(bases[finite]), np.max(bases[finite]), 10000)
        fitted = _polynomial(x_fit, *popt)

        plt.plot(bases, heights_opt, 'o', color='blue', markersize=3, label='a')
        plt.plot(x_fit, fitted, color='red', linewidth=1.0, linestyle='-')
        plt.xlabel('baseline [V]', fontsize=16)
        plt.ylabel(f'{ValueKey} [V]', fontsize=16)
        if Title is not None:
            plt.title(f"{Title} (fit)")
        plt.grid()
        plt.show()
        plt.cla()

        plt.plot(bases, data[ResultKey], 'o', color='tab:blue', markersize=0.7, label='a')
        plt.xlabel('baseline [V]', fontsize=16)
        plt.ylabel(f'{ResultKey} [V]', fontsize=16)
        if Title is not None:
            plt.title(f"{Title} (calibrated)")
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


# NoiseSPE が何であるかの明示（項目4）:
#   NoiseAnalysis() が modelnoise.txt に書き出しているのは片側 **ASD**
#   （振幅密度, unit/sqrt(Hz)。ここでは eta 換算後なので pA/sqrt(Hz)）であり、
#   FFT振幅そのものでも PSD でもない。したがって
#       PSD(f) = NoiseSPE(f) ** 2
#   Getpara の modelnoise.txt も |FFT|/sqrt(df)*eta という ASD なので、
#   「NoiseSPE で割る」= ASD 重み付け = 1/sqrt(PSD) 重み付けであり、
#   理論的な最適フィルタ S*(f)/PSD(f) とは重みが異なる。
LEGACY_FFT_METHOD = "Legacy (fft/ifft)"
CURRENT_METHOD = "Current (rfft/irfft + Bessel)"
CURRENT_NO_TEMPLATE_BESSEL_METHOD = "Current (rfft/irfft, no template Bessel)"
PSD_OPTIMAL_METHOD = "PSD-optimal (S*/PSD, normalized)"

OPTIMAL_FILTER_METHODS = (
    CURRENT_METHOD,
    CURRENT_NO_TEMPLATE_BESSEL_METHOD,
    PSD_OPTIMAL_METHOD,
    LEGACY_FFT_METHOD,
)


def _full_spectrum_sum(one_sided, sample):
    """片側rfft配列の値を、両側FFTでの総和に直して返す。

    ``sum_n a[n]*b[n] = (1/N) * sum_k^{full} A[k] conj(B[k])`` を片側配列だけで
    評価するために使う。内側のビンは正負の周波数ぶん2回数える。
    """
    one_sided = np.asarray(one_sided, dtype=float)
    if sample % 2 == 0:
        return one_sided[0] + 2.0 * np.sum(one_sided[1:-1]) + one_sided[-1]
    return one_sided[0] + 2.0 * np.sum(one_sided[1:])


def MatchedFilterNormalization(S, psd, sample, band=None):
    """``(1/N) * Σ_full |S|^2 / PSD`` を返す（振幅推定量の規格化因子）。"""
    S = np.asarray(S)
    psd = np.asarray(psd, dtype=float)
    weight = np.zeros(len(S), dtype=float)
    valid = psd > 0
    if band is not None:
        valid &= band
    weight[valid] = np.abs(S[valid]) ** 2 / psd[valid]
    return _full_spectrum_sum(weight, sample) / sample


def OptimalFilterTemplate(
    NoiseSPE,
    AveragePulse,
    config,
    method=CURRENT_METHOD,
    plot=True,
    AmplitudeScale=None,
):
    """最適フィルタのテンプレート（時間領域）を返す。

    どの方式でも ``np.sum(pulse * filt)`` で使う前提。実数波形どうしの
    内積は周波数領域では ``(1/N) Σ_full D(f) conj(H(f))`` になるので、
    ``filt = irfft(W)`` としたとき推定量は ``(1/N) Σ_full D(f) conj(W(f))``
    となる。つまり ``W = S/PSD`` と置けば推定量は自動的に
    ``Σ S*(f) D(f) / PSD(f)`` になり、複素共役は時間領域の相関側が担う
    （テンプレート側で conj を取ると符号付き虚部が二重に反転する）。

    method:
      ``"Legacy (fft/ifft)"``
          旧版と同じ全周波数FFT/ifft。重みは 1/ASD。
      ``"Current (rfft/irfft + Bessel)"``
          rfft/irfft。重みは 1/ASD。テンプレートに追加でBesselを掛ける
          （各パルスにも掛かるので帯域制限は二重になる）。
      ``"Current (rfft/irfft, no template Bessel)"``
          上と同じだがテンプレート側のBesselを外す。Getparaと同じ
          「パルスにだけBessel」になる。
      ``"PSD-optimal (S*/PSD, normalized)"``
          理論的な最適フィルタ ``S*(f)/PSD(f)``。``Σ|S|^2/PSD`` で規格化
          するので、平均パルスと同じ波形に対して推定量は 1 になる。
          ``AmplitudeScale`` を渡すとその値（通常は平均パルスの波高）を
          掛けて ``Peak`` と同じスケールに揃える。

    規格化されるのは PSD-optimal のみ。他の方式は既存の ``PeakOpt`` の
    スケールを保つため従来どおり無規格化のまま。
    """
    rate = config["Readout"]["Rate"]
    sample = config["Readout"]["Sample"]
    cf = config["Analysis"]["CutoffFrequency"]

    if method == LEGACY_FFT_METHOD:
        # 旧版と同じ全周波数FFT・旧周波数軸・ifftの組み合わせ。
        # 現行のmodelnoise.txtは片側長の場合があるため、除算できる全長へ
        # 補間する（入力ファイル形式自体は変更しない）。
        fq = np.arange(0, rate, rate / sample)
        F = scipy.fftpack.fft(AveragePulse)
        asd = _resize_spectrum(NoiseSPE, sample)
        F_filtered = np.copy(F)
        F_filtered[fq > cf] = 0
        filt_time = scipy.fftpack.ifft(F_filtered / asd).real
    elif method in (CURRENT_METHOD, CURRENT_NO_TEMPLATE_BESSEL_METHOD):
        # 現行方式：実信号用の片側FFTとirfftを使用。重みは 1/ASD。
        fq = np.fft.rfftfreq(sample, d=1 / rate)
        F = np.fft.rfft(AveragePulse)
        asd = _resize_spectrum(NoiseSPE, len(F))
        F_filtered = np.copy(F)
        F_filtered[fq > cf] = 0
        filt_time = np.fft.irfft(F_filtered / asd, n=sample).real
        if method == CURRENT_METHOD:
            # 各パルスにもBesselが掛かるので、これは二重の帯域制限になる。
            filt_time = Bessel(filt_time, rate, cf)
    elif method == PSD_OPTIMAL_METHOD:
        fq = np.fft.rfftfreq(sample, d=1 / rate)
        S = np.fft.rfft(AveragePulse)
        asd = _resize_spectrum(NoiseSPE, len(S))
        psd = np.asarray(asd, dtype=float) ** 2  # ASD -> PSD

        band = (fq <= cf) & (psd > 0)
        weight = np.zeros(len(S), dtype=complex)
        weight[band] = S[band] / psd[band]

        norm = MatchedFilterNormalization(S, psd, sample, band=band)
        if not np.isfinite(norm) or norm <= 0:
            raise ValueError(
                "OptimalFilterTemplate: Σ|S|^2/PSD が0以下です。"
                "ノイズモデルと平均パルスの整合を確認してください。"
            )
        # ここで割ると sum(AveragePulse * filt) == 1 になる。
        weight /= norm
        filt_time = np.fft.irfft(weight, n=sample)
        if AmplitudeScale is not None:
            # PeakやPeakOptと同じ単位で読めるよう、波高スケールへ戻す。
            filt_time = filt_time * float(AmplitudeScale)
    else:
        raise ValueError(f"Unknown optimal filter method: {method}")

    if plot:
        time = np.arange(sample) / rate
        plt.plot(time, filt_time)
        plt.title(f"Optimal Filter (time domain)\n{method}")
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
