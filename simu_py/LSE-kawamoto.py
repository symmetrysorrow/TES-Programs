import json
import time as pytime

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import minimize

import general as gn
import pulse_model


sim_path = "d:/kamoto"
exp_path = "G:/TSURUTA/20230616_post/room1-ch2-3_180mK_570uA_100kHz_g10"

eta = 101
amp = 100


def extract_features(pulse, dt=1.0, presample=500):
    y = pulse - np.mean(pulse[:presample])

    h_max = np.max(y)
    if h_max <= 0:
        return 0.0, 0.0, 1e-9

    idx_max = np.argmax(y)

    try:
        decay_part = y[idx_max:]

        start_candidates = np.where(decay_part <= h_max * 0.8)[0]
        end_candidates = np.where(decay_part <= h_max * 0.2)[0]

        if len(start_candidates) == 0 or len(end_candidates) == 0:
            tau = 1e-9
        else:
            start_idx = start_candidates[0]
            end_idx = end_candidates[0]

            fit_y = decay_part[start_idx:end_idx]
            fit_t = np.arange(len(fit_y)) * dt

            fit_y = fit_y[fit_y > 0]
            fit_t = fit_t[:len(fit_y)]

            if len(fit_y) > 5:
                slope, _ = np.polyfit(fit_t, np.log(fit_y), 1)
                tau = -1.0 / slope if slope < 0 else 1e-9
            else:
                tau = 1e-9

    except Exception:
        tau = 1e-9

    rise_10_candidates = np.where(y[:idx_max + 1] >= h_max * 0.1)[0]
    rise_90_candidates = np.where(y[:idx_max + 1] >= h_max * 0.9)[0]

    if len(rise_10_candidates) == 0 or len(rise_90_candidates) == 0:
        rise_time = 0.0
    else:
        rise_10 = rise_10_candidates[0]
        rise_90 = rise_90_candidates[0]
        rise_time = (rise_90 - rise_10) * dt

    return h_max, rise_time, tau


def find_decay_point_on_curve(pulse, time_axis, percent, presample):
    """
    ピーク後に波高のpercentまで減衰した点を、
    実際のプロット波形上の座標として返す。

    return:
        {
            "t": プロット上の時刻,
            "y": プロット上の振幅,
            "idx": 補間後のindex,
        }
        見つからなければ None
    """

    pulse = np.asarray(pulse)
    time_axis = np.asarray(time_axis)

    baseline = np.mean(pulse[:presample])
    y = pulse - baseline

    h_max = np.max(y)
    if h_max <= 0:
        return None

    idx_peak = int(np.argmax(y))
    target = h_max * percent

    for i in range(idx_peak + 1, len(y)):
        if y[i] <= target:
            y1 = y[i - 1]
            y2 = y[i]

            if y1 == y2:
                idx_float = float(i)
            else:
                idx_float = (i - 1) + (target - y1) / (y2 - y1)

            t_float = np.interp(
                idx_float,
                np.arange(len(time_axis)),
                time_axis,
            )

            y_plot = baseline + target

            return {
                "t": t_float,
                "y": y_plot,
                "idx": idx_float,
            }

    return None


def extract_decay_points_on_curve(pulse, time_axis, presample):
    p90 = find_decay_point_on_curve(
        pulse=pulse,
        time_axis=time_axis,
        percent=0.9,
        presample=presample,
    )

    p50 = find_decay_point_on_curve(
        pulse=pulse,
        time_axis=time_axis,
        percent=0.5,
        presample=presample,
    )

    return p90, p50


def MakePulse(para, Presample):
    ch0, ch1 = pulse_model.model(para)
    pulses = []

    for pulse in ch0:
        zeros = np.zeros(int(Presample))
        pulse = np.concatenate([zeros, pulse])[:int(para["samples"])]
        pulses.append(pulse)

    return pulses


def normalize_pulse(pulse, presample):
    y = pulse.copy()
    y -= np.mean(y[:presample])
    return y


def build_wave_weights(pulse, presample):
    y = normalize_pulse(pulse, presample)
    weights = np.full(len(y), 0.2)

    h_max = np.max(y)
    if h_max <= 0:
        return weights

    idx_peak = np.argmax(y)

    above_10 = np.where(y[:idx_peak + 1] >= h_max * 0.1)[0]
    above_90 = np.where(y[:idx_peak + 1] >= h_max * 0.9)[0]

    rise_start = int(above_10[0]) if len(above_10) else max(0, idx_peak - 5)
    rise_end = int(above_90[0]) if len(above_90) else idx_peak
    rise_width = max(1, rise_end - rise_start)

    peak_core_left = max(0, idx_peak - rise_width)
    peak_core_right = min(len(y), idx_peak + rise_width + 1)

    peak_shoulder_left = max(0, idx_peak - 3 * rise_width)
    peak_shoulder_right = min(len(y), idx_peak + 3 * rise_width + 1)

    decay_tail_right = min(len(y), idx_peak + 8 * rise_width + 1)

    weights[peak_shoulder_left:peak_shoulder_right] = np.maximum(
        weights[peak_shoulder_left:peak_shoulder_right],
        6.0,
    )

    weights[peak_core_left:peak_core_right] = np.maximum(
        weights[peak_core_left:peak_core_right],
        20.0,
    )

    weights[idx_peak] = 30.0

    weights[idx_peak:decay_tail_right] = np.maximum(
        weights[idx_peak:decay_tail_right],
        1.5,
    )

    return weights


def weighted_wave_error(sim_pulse, exp_pulse, weights, presample):
    sim = normalize_pulse(sim_pulse, presample)
    exp = normalize_pulse(exp_pulse, presample)

    scale = max(np.max(np.abs(exp)), 1e-12)

    return np.sum(weights * ((sim - exp) / scale) ** 2) / np.sum(weights)


def plot_pulses_with_points(
    time_axis,
    pulse_high,
    pulse_low,
    pulse_sim_high,
    pulse_sim_low,
    Presample,
    title,
):
    exp_p90, exp_p50 = extract_decay_points_on_curve(
        pulse=pulse_high,
        time_axis=time_axis,
        presample=Presample,
    )

    sim_p90, sim_p50 = extract_decay_points_on_curve(
        pulse=pulse_sim_high,
        time_axis=time_axis,
        presample=Presample,
    )

    plt.figure(figsize=(12, 7))

    plt.plot(
        time_axis,
        pulse_high,
        label="Experimental High Pulse",
        color="gray",
        linestyle="--",
    )

    plt.plot(
        time_axis,
        pulse_low,
        label="Experimental Low Pulse",
        color="lightgray",
        linestyle="--",
    )

    plt.plot(
        time_axis,
        pulse_sim_high,
        label="Simulated High Pulse",
    )

    plt.plot(
        time_axis,
        pulse_sim_low,
        label="Simulated Low Pulse",
    )

    if exp_p90 is not None:
        plt.scatter(
            exp_p90["t"],
            exp_p90["y"],
            s=100,
            marker="o",
            label="Exp High 90%",
            zorder=20,
        )

    if exp_p50 is not None:
        plt.scatter(
            exp_p50["t"],
            exp_p50["y"],
            s=100,
            marker="o",
            label="Exp High 50%",
            zorder=20,
        )

    if sim_p90 is not None:
        plt.scatter(
            sim_p90["t"],
            sim_p90["y"],
            s=120,
            marker="x",
            label="Sim High 90%",
            zorder=20,
        )

    if sim_p50 is not None:
        plt.scatter(
            sim_p50["t"],
            sim_p50["y"],
            s=120,
            marker="x",
            label="Sim High 50%",
            zorder=20,
        )

    plt.xlabel("Time (s)")
    plt.ylabel("Amplitude")
    plt.title(title)
    plt.grid(True)
    plt.legend()
    plt.show()

    return exp_p90, exp_p50, sim_p90, sim_p50


with open(f"{sim_path}/input.json", "r") as f:
    para = json.load(f)


selected_keys = gn.LoadTxt(f"{exp_path}/SelectedKeys_fromScatter.txt")
para_exp = gn.LoadJson(f"{exp_path}/PulseConfig.json")

Presample = int(para_exp["Readout"]["PreSample"])
dt = para_exp["Readout"]["Sample"] / para_exp["Readout"]["Rate"]

time_axis = gn.GetTime(
    para_exp["Readout"]["Rate"],
    para_exp["Readout"]["Sample"],
)


peak_values = []
keys_list = []
total_keys = len(selected_keys)

for i, key in enumerate(selected_keys, start=1):
    pulse = gn.LoadBin(f"{exp_path}/CH0_pulse/rawdata/CH0_{int(key)}.dat")

    pulse = pulse.copy()
    pulse -= np.mean(pulse[0:Presample])

    pulse = gn.Bessel(
        pulse,
        para_exp["Readout"]["Rate"],
        para_exp["Analysis"]["CutoffFrequency"],
    )

    peak_values.append(np.max(pulse))
    keys_list.append(key)

    if i == 1 or i == total_keys or i % max(1, total_keys // 10) == 0:
        print(f"[Load] {i}/{total_keys} pulses processed")


peak_values = np.array(peak_values)

max_idx = np.argmax(peak_values)
min_idx = np.argmin(peak_values)

key_high = keys_list[max_idx]
key_low = keys_list[min_idx]


pulse_high = (
    gn.LoadBin(f"{exp_path}/CH0_pulse/rawdata/CH0_{int(key_high)}.dat")
    * eta
    / amp
    / 1e6
)
pulse_high -= np.mean(pulse_high[0:Presample])

pulse_low = (
    gn.LoadBin(f"{exp_path}/CH0_pulse/rawdata/CH0_{int(key_low)}.dat")
    * eta
    / amp
    / 1e6
)
pulse_low -= np.mean(pulse_low[0:Presample])


exp_h_high, exp_tr_high, exp_tau = extract_features(
    pulse_high,
    dt=dt,
    presample=Presample,
)

exp_h_low, exp_tr_low, _ = extract_features(
    pulse_low,
    dt=dt,
    presample=Presample,
)


exp_p90_high, exp_p50_high = extract_decay_points_on_curve(
    pulse=pulse_high,
    time_axis=time_axis,
    presample=Presample,
)

print(f"Highest Pulse: Height = {exp_h_high}, Rise Time = {exp_tr_high}, Tau = {exp_tau}")
print(f"Lowest Pulse : Height = {exp_h_low}, Rise Time = {exp_tr_low}")

if exp_p90_high is not None:
    print(f"Experimental High Pulse t90 = {exp_p90_high['t']}")
else:
    print("Experimental High Pulse t90 = None")

if exp_p50_high is not None:
    print(f"Experimental High Pulse t50 = {exp_p50_high['t']}")
else:
    print("Experimental High Pulse t50 = None")


wave_weights_high = build_wave_weights(pulse_high, Presample)
wave_weights_low = build_wave_weights(pulse_low, Presample)


pulses = MakePulse(para, Presample)

pulse_sim_high = pulses[0]
pulse_sim_low = pulses[-1]

init_h_high, init_tr_high, init_tau_high = extract_features(
    pulse_sim_high,
    dt=dt,
    presample=Presample,
)

init_h_low, init_tr_low, init_tau_low = extract_features(
    pulse_sim_low,
    dt=dt,
    presample=Presample,
)

init_p90_high, init_p50_high = extract_decay_points_on_curve(
    pulse=pulse_sim_high,
    time_axis=time_axis,
    presample=Presample,
)

print(
    f"Initial Simulation High Pulse: "
    f"Height = {init_h_high}, Rise Time = {init_tr_high}, Tau = {init_tau_high}"
)

print(
    f"Initial Simulation Low Pulse : "
    f"Height = {init_h_low}, Rise Time = {init_tr_low}, Tau = {init_tau_low}"
)

if init_p90_high is not None:
    print(f"Initial Simulation High Pulse t90 = {init_p90_high['t']}")
else:
    print("Initial Simulation High Pulse t90 = None")

if init_p50_high is not None:
    print(f"Initial Simulation High Pulse t50 = {init_p50_high['t']}")
else:
    print("Initial Simulation High Pulse t50 = None")


plot_pulses_with_points(
    time_axis=time_axis,
    pulse_high=pulse_high,
    pulse_low=pulse_low,
    pulse_sim_high=pulse_sim_high,
    pulse_sim_low=pulse_sim_low,
    Presample=Presample,
    title="Initial Pulse Comparison",
)


progress = {
    "eval_count": 0,
    "callback_count": 0,
    "start_time": pytime.time(),
    "last_error": None,
}


def err_func(params):
    progress["eval_count"] += 1

    if any(p <= 0 for p in params):
        return 1e18

    para["C_abs"] = params[0]
    para["C_tes"] = params[1]
    para["G_abs-abs"] = params[2]
    para["G_abs-tes"] = params[3]
    para["G_tes-bath"] = params[4]

    try:
        pulses = MakePulse(para, Presample)

        pulse_sim_high = pulses[0]
        pulse_sim_low = pulses[-1]

        sim_h_low, sim_tr_low, _ = extract_features(
            pulse_sim_low,
            dt=dt,
            presample=Presample,
        )

        sim_h_high, sim_tr_high, sim_tau = extract_features(
            pulse_sim_high,
            dt=dt,
            presample=Presample,
        )

        sim_p90_high, sim_p50_high = extract_decay_points_on_curve(
            pulse=pulse_sim_high,
            time_axis=time_axis,
            presample=Presample,
        )

    except Exception:
        return 1e18

    if (
        exp_p90_high is None
        or exp_p50_high is None
        or sim_p90_high is None
        or sim_p50_high is None
    ):
        return 1e18

    exp_t90_high = exp_p90_high["t"]
    exp_t50_high = exp_p50_high["t"]

    sim_t90_high = sim_p90_high["t"]
    sim_t50_high = sim_p50_high["t"]

    e_wave_high = weighted_wave_error(
        pulse_sim_high,
        pulse_high,
        wave_weights_high,
        Presample,
    )

    e_wave_low = weighted_wave_error(
        pulse_sim_low,
        pulse_low,
        wave_weights_low,
        Presample,
    )

    e_h_high = ((sim_h_high - exp_h_high) / exp_h_high) ** 2
    e_tr_high = ((sim_tr_high - exp_tr_high) / exp_tr_high) ** 2

    e_h_low = ((sim_h_low - exp_h_low) / exp_h_low) ** 2
    e_tr_low = ((sim_tr_low - exp_tr_low) / exp_tr_low) ** 2

    e_tau = ((sim_tau - exp_tau) / exp_tau) ** 2

    e_t90_high = ((sim_t90_high - exp_t90_high) / exp_t90_high) ** 2
    e_t50_high = ((sim_t50_high - exp_t50_high) / exp_t50_high) ** 2

    wave_errors = np.array([
        e_wave_high,
        e_wave_low,
    ])

    wave_error_weights = np.array([
        30.0,
        20.0,
    ])

    feature_errors = np.array([
        e_h_high,
        e_tr_high,
        e_h_low,
        e_tr_low,
        e_tau,
        e_t90_high,
        e_t50_high,
    ])

    feature_error_weights = np.array([
        8.0,
        12.0,
        8.0,
        6.0,
        6.0,
        10.0,
        10.0,
    ])

    total_error = (
        np.sum(wave_error_weights * wave_errors)
        + np.sum(feature_error_weights * feature_errors)
    )

    progress["last_error"] = total_error

    if progress["eval_count"] == 1 or progress["eval_count"] % 20 == 0:
        elapsed = pytime.time() - progress["start_time"]

        print(
            f"[Optimize] eval={progress['eval_count']}, "
            f"error={total_error:.6e}, "
            f"wave_high={e_wave_high:.3e}, wave_low={e_wave_low:.3e}, "
            f"t90_exp={exp_t90_high:.6e}, t90_sim={sim_t90_high:.6e}, "
            f"t50_exp={exp_t50_high:.6e}, t50_sim={sim_t50_high:.6e}, "
            f"elapsed={elapsed:.1f}s"
        )

    return total_error


def optimization_callback(xk):
    progress["callback_count"] += 1

    elapsed = pytime.time() - progress["start_time"]

    params_str = ", ".join(f"{v:.3e}" for v in xk)

    error_str = (
        f"{progress['last_error']:.6e}"
        if progress["last_error"] is not None
        else "N/A"
    )

    print(
        f"[Optimize] iter={progress['callback_count']}, "
        f"last_error={error_str}, elapsed={elapsed:.1f}s"
    )
    print(f"           params=[{params_str}]")


initial_params = [
    7.9e-10,
    7.9e-12,
    1.5e-07,
    8.2e-09,
    1.68e-08,
]


print("[Optimize] Starting Nelder-Mead optimization")

result = minimize(
    err_func,
    initial_params,
    method="Nelder-Mead",
    callback=optimization_callback,
    options={
        "maxiter": 500,
        "disp": True,
        "xatol": 1e-8,
        "fatol": 1e-8,
    },
)


print("Optimization Result:")
print(result.x)


para["C_abs"] = result.x[0]
para["C_tes"] = result.x[1]
para["G_abs-abs"] = result.x[2]
para["G_abs-tes"] = result.x[3]
para["G_tes-bath"] = result.x[4]


with open(f"{sim_path}/final_params.json", "w") as f:
    json.dump(para, f, indent=4)


pulses = MakePulse(para, Presample)

final_high = pulses[0]
final_low = pulses[-1]

final_h_high, final_tr_high, final_tau_high = extract_features(
    final_high,
    dt=dt,
    presample=Presample,
)

final_h_low, final_tr_low, final_tau_low = extract_features(
    final_low,
    dt=dt,
    presample=Presample,
)

final_p90_high, final_p50_high = extract_decay_points_on_curve(
    pulse=final_high,
    time_axis=time_axis,
    presample=Presample,
)


print(
    f"High Pulse Simulated: "
    f"Height = {final_h_high}, Rise Time = {final_tr_high}, Tau = {final_tau_high}"
)

print(
    f"Low Pulse Simulated : "
    f"Height = {final_h_low}, Rise Time = {final_tr_low}, Tau = {final_tau_low}"
)

if final_p90_high is not None:
    print(f"Final Simulation High Pulse t90 = {final_p90_high['t']}")
else:
    print("Final Simulation High Pulse t90 = None")

if final_p50_high is not None:
    print(f"Final Simulation High Pulse t50 = {final_p50_high['t']}")
else:
    print("Final Simulation High Pulse t50 = None")


plot_pulses_with_points(
    time_axis=time_axis,
    pulse_high=pulse_high,
    pulse_low=pulse_low,
    pulse_sim_high=final_high,
    pulse_sim_low=final_low,
    Presample=Presample,
    title="Final Pulse Comparison",
)