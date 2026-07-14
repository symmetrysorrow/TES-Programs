import numpy as np
import json
import pandas as pd
import os
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from scipy.interpolate import UnivariateSpline
import natsort
import glob
import re
import tqdm
import getpara as gp


def save_readpulse_debug(Data_path, pulse_raw, pulse_filt, para, path, reason, peak_index=None, rise_10=None, rise_90=None):
    debug_dir = f"{Data_path}/debug_readpulse"
    os.makedirs(debug_dir, exist_ok=True)

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
    plt.savefig(f"{debug_dir}/{safe_path}.png", dpi=200)
    plt.close()


def ReadOutput(FilePath, key):
    df = pd.read_csv(FilePath)
    return df[key].to_numpy()


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

    sum_height = CH0_df["height"].to_numpy(dtype=float) + CH1_df["height"].to_numpy(dtype=float)
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


def remove_outliers(data, percentile=0.1):
    lower_bound = np.percentile(data, percentile)
    upper_bound = np.percentile(data, 100 - percentile)
    return data[(data >= lower_bound) & (data <= upper_bound)]


def ReadPulse(Data_path, pulse, target="Pulse_ms", path="", debug_plot=False):
    with open(f"{Data_path}/input.json") as f:
        para = json.load(f)

    if np.mean(pulse) <= 0:
        pulse *= -1

    pulse_filt = gp.BesselFilter(pulse, para["rate"], para["cutoff"]) if target != "Pulse_ms" else pulse

    if not np.all(np.isfinite(pulse_filt)):
        print(f"ReadPulse failed: non-finite values in pulse ({path})")
        save_readpulse_debug(Data_path, pulse, pulse_filt, para, path, "non-finite values")
        return [np.nan, np.nan, np.nan, np.nan]

    peak = np.max(pulse_filt)
    peak_index = int(np.argmax(pulse_filt))

    if peak_index < 10:
        print(f"ReadPulse failed: peak too early ({path})")
        save_readpulse_debug(Data_path, pulse, pulse_filt, para, path, "peak too early", peak_index=peak_index)
        return [np.nan, np.nan, np.nan, np.nan]

    l = max(0, peak_index - 10)
    r = min(len(pulse_filt), peak_index + 90)
    peak_av = np.mean(pulse_filt[l:r])

    rise_90 = None
    for i in reversed(range(0, peak_index)):
        if pulse_filt[i] <= peak * 0.9:
            rise_90 = i
            break
    if rise_90 is None:
        print(f"ReadPulse failed: rise_90 not found ({path})")
        save_readpulse_debug(Data_path, pulse, pulse_filt, para, path, "rise_90 not found", peak_index=peak_index)
        return [np.nan, np.nan, np.nan, np.nan]

    rise_10 = None
    for j in reversed(range(0, rise_90)):
        if pulse_filt[j] <= peak * 0.1:
            rise_10 = j
            break
    if rise_10 is None:
        print(f"ReadPulse failed: rise_10 not found ({path})")
        save_readpulse_debug(Data_path, pulse, pulse_filt, para, path, "rise_10 not found", peak_index=peak_index, rise_90=rise_90)
        return [np.nan, np.nan, np.nan, np.nan]

    rise = (rise_90 - rise_10) / para["rate"]

    st_l = para["SettlingTime"] - 10
    st_r = para["SettlingTime"] + 90
    if st_l < 0 or st_r > len(pulse_filt) or st_l >= st_r:
        print(f"ReadPulse failed: invalid settling window ({path})")
        save_readpulse_debug(Data_path, pulse, pulse_filt, para, path, "invalid settling window", peak_index=peak_index, rise_10=rise_10, rise_90=rise_90)
        return [np.nan, np.nan, np.nan, np.nan]

    ST_window = pulse_filt[st_l:st_r]
    if len(ST_window) == 0 or not np.all(np.isfinite(ST_window)):
        print(f"ReadPulse failed: invalid settling samples ({path})")
        save_readpulse_debug(Data_path, pulse, pulse_filt, para, path, "invalid settling samples", peak_index=peak_index, rise_10=rise_10, rise_90=rise_90)
        return [np.nan, np.nan, np.nan, np.nan]

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

    return [peak_av, peak_index, rise, ST_height]


def MakeOutput(Data_path, target):
    with open(f"{Data_path}/input.json") as f:
        para = json.load(f)

    if target != "Pulse_ms":
        print("BesselFilter")

    for i, posi in enumerate(para["position"]):
        print(f"{i+1}/{len(para['position'])}")
        for ch in [0, 1]:
            print(f"Ch:{ch}")
            results = []
            pulse_numbers = []
            pulse_pathes = natsort.natsorted(glob.glob(f'{Data_path}/{para["E"]}keV_{posi}/{target}/CH{ch}/CH{ch}_*.dat'))

            for path in tqdm.tqdm(pulse_pathes):
                match = re.search(fr'CH{ch}_(\d+).dat', path)
                pulse_numbers.append(match.group(1))
                pulse = np.loadtxt(path)
                result = ReadPulse(Data_path, pulse, path=path, target=target)
                results.append(result)

            columns = ["height", "peak_index", "rise", "ST_Height"]
            df = pd.DataFrame(results, columns=columns, index=pulse_numbers)
            df.index.name = "id"
            df.to_csv(f'{Data_path}/{para["E"]}keV_{posi}/{target}/output_TES{ch}.csv')


def gaussian(x, amp, mean, stddev):
    return amp * np.exp(-((x - mean) ** 2) / (2 * stddev ** 2))


def find_nearest_values(ratios, Fit_para):
    B_right = Fit_para[:, 1]
    B_left = Fit_para[:, 0]
    result = []
    for a in ratios:
        idx = np.abs(B_right - a).argmin()
        result.append(B_left[idx])
    return np.array(result)


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


def MakeHistgram(data, posi, HistColor=None):
    data = np.asarray(data, dtype=float)
    data = data[np.isfinite(data)]
    if len(data) == 0:
        return np.nan, np.nan

    bin_num = optimal_bin_count(data)
    hist, bin_edges = np.histogram(data, bins=bin_num, density=False)
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    initial_guess = [np.max(hist), np.mean(data), np.std(data)]

    if HistColor is not None:
        plt.hist(data, bins=bin_num, density=False, label=f"abs-{posi}", color=HistColor)
    else:
        plt.hist(data, bins=bin_num, density=False, label=f"abs-{posi}")

    popt, _ = curve_fit(gaussian, bin_centers, hist, p0=initial_guess, maxfev=10000000)
    _, mean_fit, stddev_fit = popt
    fwhm = 2 * stddev_fit * np.sqrt(2 * np.log(2))

    reference_width = 2 * np.std(data) * np.sqrt(2 * np.log(2))
    sample_mean = np.mean(data)
    robust_sigma = (np.percentile(data, 84) - np.percentile(data, 16)) / 2
    robust_fwhm = 2 * robust_sigma * np.sqrt(2 * np.log(2))
    robust_reso = robust_fwhm / sample_mean if sample_mean > 0 else np.nan
    fit_reso = fwhm / mean_fit if mean_fit > 0 else np.nan

    if np.isfinite(reference_width) and reference_width > 0 and fwhm < 0.6 * reference_width:
        retry_bins = min(max(12, bin_num // 2), 14)
        if retry_bins != bin_num:
            hist, bin_edges = np.histogram(data, bins=retry_bins, density=False)
            bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
            initial_guess = [np.max(hist), np.mean(data), np.std(data)]
            popt, _ = curve_fit(gaussian, bin_centers, hist, p0=initial_guess, maxfev=10000000)
            _, mean_fit, stddev_fit = popt
            fwhm = 2 * stddev_fit * np.sqrt(2 * np.log(2))
            fit_reso = fwhm / mean_fit if mean_fit > 0 else np.nan

    # If the Gaussian fit becomes much broader than a robust percentile-based
    # estimate, fall back to the robust width so one skewed histogram binning
    # does not dominate the reported resolution.
    if np.isfinite(robust_reso) and np.isfinite(fit_reso) and fit_reso > 1.5 * robust_reso:
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


def Resos(Data_path, target, show):
    fit_para = np.loadtxt(f"{Data_path}/ratios.txt", delimiter=',')

    with open(f"{Data_path}/input.json", "r") as f:
        para = json.load(f)

    CH0_files = []
    CH1_files = []
    for posi in para["position"]:
        CH0_files.append(f"{Data_path}/{para['E']}keV_{posi}/{target}/output_TES0.csv")
        CH1_files.append(f"{Data_path}/{para['E']}keV_{posi}/{target}/output_TES1.csv")

    ratio_base = fit_para[:, 1]
    posi_base = fit_para[:, 0]
    sorted_indices = np.argsort(ratio_base)
    ratio_base = ratio_base[sorted_indices]
    posi_base = posi_base[sorted_indices]
    spline = UnivariateSpline(ratio_base, posi_base, k=3, s=0)

    fwhms = []
    for CH0, CH1, posi in zip(CH0_files, CH1_files, para["position"]):
        CH0_df = pd.read_csv(CH0)
        CH1_df = pd.read_csv(CH1)
        CH0_df, CH1_df, quality_mask = remove_outlier_events(CH0_df, CH1_df)

        CH0_heights = CH0_df["height"].to_numpy(dtype=float)
        CH1_heights = CH1_df["height"].to_numpy(dtype=float)
        mask = np.isfinite(CH0_heights) & np.isfinite(CH1_heights) & (CH1_heights != 0)
        ratios = CH0_heights[mask] / CH1_heights[mask]

        positions = np.asarray(spline(ratios), dtype=float)
        positions = positions[np.isfinite(positions)]

        fwhm, _ = MakeHistgram(positions, posi)
        fwhms.append(fwhm)
        removed = int((~quality_mask).sum())
        if removed > 0:
            print(f"Position {posi}: removed {removed} outlier events from CH0/CH1 pair")

    fwhms = np.array(fwhms)
    np.savetxt(f"{Data_path}/fwhms_{target}.txt", fwhms)

    plt.xlabel("Position", fontsize=15)
    plt.ylabel("Count", fontsize=15)
    plt.title(f"{para['E']}keV")
    plt.tight_layout()
    plt.legend()
    plt.savefig(f"{Data_path}/position_histgram_{target}_{para['E']}.png")
    if show:
        plt.show()
        plt.plot(para["position"], fwhms)
        plt.show()
    else:
        plt.clf()

    colors = generate_symmetric_colors(len(fwhms))

    ene_reso_sum = []
    for CH0, CH1, posi, color in zip(CH0_files, CH1_files, para["position"], colors):
        CH0_df = pd.read_csv(CH0)
        CH1_df = pd.read_csv(CH1)
        CH0_df, CH1_df, quality_mask = remove_outlier_events(CH0_df, CH1_df)
        CH0_heights = CH0_df["height"].to_numpy(dtype=float)
        CH1_heights = CH1_df["height"].to_numpy(dtype=float)
        _, reso = MakeHistgram(CH0_heights + CH1_heights, posi, color)
        ene_reso_sum.append(reso * para["E"])
        removed = int((~quality_mask).sum())
        if removed > 0:
            print(f"Position {posi}: removed {removed} outlier events for Sum")

    plt.xlabel("Current[A]", fontsize=15)
    plt.ylabel("Count", fontsize=15)
    plt.tight_layout()
    plt.legend(labelspacing=0, fontsize=8, markerscale=0.5)
    plt.savefig(f"{Data_path}/energy_sum_histgram_{target}_{para['E']}.png")
    if show:
        plt.show()
        plt.plot(para["position"], ene_reso_sum)
        plt.show()
    else:
        plt.clf()

    ene_reso_max = []
    for CH0, CH1, posi, color in zip(CH0_files, CH1_files, para["position"], colors):
        CH0_df = pd.read_csv(CH0)
        CH1_df = pd.read_csv(CH1)
        CH0_df, CH1_df, quality_mask = remove_outlier_events(CH0_df, CH1_df)
        CH0_heights = CH0_df["height"].to_numpy(dtype=float)
        CH1_heights = CH1_df["height"].to_numpy(dtype=float)
        _, reso = MakeHistgram(np.maximum(CH0_heights, CH1_heights), posi, color)
        ene_reso_max.append(reso * para["E"])
        removed = int((~quality_mask).sum())
        if removed > 0:
            print(f"Position {posi}: removed {removed} outlier events for Max")

    plt.xlabel("Current[A]", fontsize=15)
    plt.ylabel("Count", fontsize=15)
    plt.tight_layout()
    plt.xlim(0, None)
    plt.legend(labelspacing=0, fontsize=8, markerscale=0.5)
    plt.savefig(f"{Data_path}/energy_max_histgram_{target}_{para['E']}.png")
    if show:
        plt.show()
        plt.plot(para["position"], ene_reso_max)
        plt.show()
    else:
        plt.clf()

    ene_reso_min = []
    for CH0, CH1, posi, color in zip(CH0_files, CH1_files, para["position"], colors):
        CH0_df = pd.read_csv(CH0)
        CH1_df = pd.read_csv(CH1)
        CH0_df, CH1_df, quality_mask = remove_outlier_events(CH0_df, CH1_df)
        CH0_heights = CH0_df["height"].to_numpy(dtype=float)
        CH1_heights = CH1_df["height"].to_numpy(dtype=float)
        _, reso = MakeHistgram(np.minimum(CH0_heights, CH1_heights), posi, color)
        ene_reso_min.append(reso * para["E"])
        removed = int((~quality_mask).sum())
        if removed > 0:
            print(f"Position {posi}: removed {removed} outlier events for Min")

    plt.xlabel("Current[A]", fontsize=15)
    plt.ylabel("Count", fontsize=15)
    plt.tight_layout()
    plt.legend(labelspacing=0, fontsize=8, markerscale=0.5)
    plt.savefig(f"{Data_path}/energy_min_histgram_{target}_{para['E']}.png")
    if show:
        plt.show()
        plt.plot(para["position"], ene_reso_min)
        plt.show()
    else:
        plt.clf()

    ene_reso_st = []
    for CH0, CH1, posi, color in zip(CH0_files, CH1_files, para["position"], colors):
        CH0_df = pd.read_csv(CH0)
        CH1_df = pd.read_csv(CH1)
        CH0_df, CH1_df, quality_mask = remove_outlier_events(CH0_df, CH1_df)
        CH0_heights = CH0_df["ST_Height"].to_numpy(dtype=float)
        CH1_heights = CH1_df["ST_Height"].to_numpy(dtype=float)
        _, reso = MakeHistgram(CH0_heights + CH1_heights, posi, color)
        ene_reso_st.append(reso * para["E"])
        removed = int((~quality_mask).sum())
        if removed > 0:
            print(f"Position {posi}: removed {removed} outlier events for ST")

    plt.xlabel("Current[A]", fontsize=15)
    plt.ylabel("Count", fontsize=15)
    plt.tight_layout()
    plt.legend(labelspacing=0, fontsize=8, markerscale=0.5)
    plt.savefig(f"{Data_path}/energy_ch_histgram_{target}_{para['E']}.png")
    if show:
        plt.show()
        plt.plot(para["position"], ene_reso_st)
        plt.show()
    else:
        plt.clf()

    data = {
        "Sum": ene_reso_sum,
        "Max": ene_reso_max,
        "Min": ene_reso_min,
        "ST": ene_reso_st,
    }
    index_label = (np.array(para["position"]) - 0.5) * para["length"] / para["n_abs"]
    df = pd.DataFrame(data, index=index_label)
    df.to_csv(f"{Data_path}/ene_resos_{target}.csv")


def try_ReadPulse(Data_path):
    pulse_path = input("Pulse path:")
    pulse = np.loadtxt(pulse_path)
    peak_av, peak_index, rise, ST_height = ReadPulse(Data_path, pulse, path=pulse_path, debug_plot=True)
    print(f"peak_av: {peak_av}, peak_index: {peak_index}, rise: {rise}, ST_height: {ST_height}")


if __name__ == "__main__":
    show = False
    Data_path = "H:\\hata\\1332_142_136_300split"

    #MakeOutput(Data_path, "Pulse_noise")
    #MakeOutput(Data_path, "Pulse_ms")
    #MakeOutput(Data_path, "Pulse_ms_noise")
    #MakeOutput(Data_path, "pulse_noise_ms_test")

    #Resos(Data_path, "Pulse_noise", show)
    Resos(Data_path, "Pulse_ms", show)
    #Resos(Data_path, "Pulse_ms_noise", show)
    #Resos(Data_path, "pulse_noise_ms_test", show)

    #try_ReadPulse(Data_path)
