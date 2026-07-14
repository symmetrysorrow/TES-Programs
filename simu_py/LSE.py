import json
import time as pytime

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import minimize

import general as gn
import pulse_model

sim_path = "h:/hata2025/200_180"
exp_path = "G:/TSURUTA/20230616_post/room1-ch2-3_180mK_570uA_100kHz_g10"
eta = 101
amp = 100
dt = 1.0


def extract_features(pulse, dt=1.0):
    y = pulse - np.mean(pulse[:500])
    h_max = np.max(y)
    if h_max <= 0:
        return 0.0, 0.0, 1e-9

    idx_max = np.argmax(y)

    try:
        decay_part = y[idx_max:]
        start_idx = np.where(decay_part <= h_max * 0.8)[0][0]
        end_idx = np.where(decay_part <= h_max * 0.2)[0][0]

        fit_y = decay_part[start_idx:end_idx]
        fit_t = np.arange(len(fit_y)) * dt

        if len(fit_y) > 5:
            slope, _ = np.polyfit(fit_t, np.log(fit_y), 1)
            tau = -1.0 / slope if slope < 0 else 1e-9
        else:
            tau = 1e-9
    except Exception:
        tau = 1e-9

    rise_20 = np.argmax(y[: idx_max + 1] >= h_max * 0.2)
    rise_80 = np.argmax(y[: idx_max + 1] >= h_max * 0.8)
    rise_time = (rise_80 - rise_20) * dt / len(pulse)
    return h_max, rise_time, tau


def MakePulse(para, Presample):
    ch0, ch1 = pulse_model.model(para)
    pulses = []
    for pulse in ch0:
        zeros = np.zeros(int(Presample))
        pulse = np.concatenate([zeros, pulse])[: int(para["samples"])]
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
    above_20 = np.where(y[: idx_peak + 1] >= h_max * 0.2)[0]
    above_80 = np.where(y[: idx_peak + 1] >= h_max * 0.8)[0]

    rise_start = int(above_20[0]) if len(above_20) else max(0, idx_peak - 5)
    rise_end = int(above_80[0]) if len(above_80) else idx_peak
    rise_width = max(1, rise_end - rise_start)

    peak_core_left = max(0, idx_peak - rise_width)
    peak_core_right = min(len(y), idx_peak + rise_width + 1)
    peak_shoulder_left = max(0, idx_peak - 3 * rise_width)
    peak_shoulder_right = min(len(y), idx_peak + 3 * rise_width + 1)
    decay_tail_right = min(len(y), idx_peak + 8 * rise_width + 1)

    weights[peak_shoulder_left:peak_shoulder_right] = np.maximum(
        weights[peak_shoulder_left:peak_shoulder_right], 6.0
    )
    weights[peak_core_left:peak_core_right] = np.maximum(
        weights[peak_core_left:peak_core_right], 20.0
    )
    weights[idx_peak] = 30.0
    weights[idx_peak:decay_tail_right] = np.maximum(weights[idx_peak:decay_tail_right], 1.5)

    return weights


def weighted_wave_error(sim_pulse, exp_pulse, weights, presample):
    sim = normalize_pulse(sim_pulse, presample)
    exp = normalize_pulse(exp_pulse, presample)
    scale = max(np.max(np.abs(exp)), 1e-12)
    return np.sum(weights * ((sim - exp) / scale) ** 2) / np.sum(weights)


with open(f"{sim_path}/input.json", "r") as f:
    para = json.load(f)

selected_keys = gn.LoadTxt(f"{exp_path}/SelectedKeys_fromScatter.txt")
para_exp = gn.LoadJson(f"{exp_path}/PulseConfig.json")
Presample = int(para_exp["Readout"]["PreSample"])
dt = para_exp["Readout"]["Sample"] / para_exp["Readout"]["Rate"]

peak_values = []
keys_list = []
total_keys = len(selected_keys)

for i, key in enumerate(selected_keys, start=1):
    pulse = gn.LoadBin(f"{exp_path}/CH0_pulse/rawdata/CH0_{int(key)}.dat")
    pulse = pulse.copy()
    pulse -= np.mean(pulse[0:Presample])
    pulse = gn.Bessel(pulse, para_exp["Readout"]["Rate"], para_exp["Analysis"]["CutoffFrequency"])

    peak_values.append(np.max(pulse))
    keys_list.append(key)

    if i == 1 or i == total_keys or i % max(1, total_keys // 10) == 0:
        print(f"[Load] {i}/{total_keys} pulses processed")

peak_values = np.array(peak_values)
max_idx = np.argmax(peak_values)
min_idx = np.argmin(peak_values)

key_high = keys_list[max_idx]
key_low = keys_list[min_idx]

pulse_high = gn.LoadBin(f"{exp_path}/CH0_pulse/rawdata/CH0_{int(key_high)}.dat") * eta / amp / 1e6
pulse_high -= np.mean(pulse_high[0:Presample])
pulse_low = gn.LoadBin(f"{exp_path}/CH0_pulse/rawdata/CH0_{int(key_low)}.dat") * eta / amp / 1e6
pulse_low -= np.mean(pulse_low[0:Presample])

exp_h_high, exp_tr_high, exp_tau = extract_features(pulse_high, dt)
exp_h_low, exp_tr_low, _ = extract_features(pulse_low, dt)
wave_weights_high = build_wave_weights(pulse_high, Presample)
wave_weights_low = build_wave_weights(pulse_low, Presample)

print(f"Highest Pulse: Height = {exp_h_high}, Rise Time = {exp_tr_high}, Tau = {exp_tau}")
print(f"Lowest Pulse : Height = {exp_h_low}, Rise Time = {exp_tr_low}")

time_axis = gn.GetTime(para_exp["Readout"]["Rate"], para_exp["Readout"]["Sample"])

pulses = MakePulse(para, Presample)
pulse_sim_high = pulses[0]
pulse_sim_low = pulses[-1]

init_h_high, init_tr_high, init_tau_high = extract_features(pulse_sim_high, dt)
init_h_low, init_tr_low, init_tau_low = extract_features(pulse_sim_low, dt)
print(f"Initial Simulation High Pulse: Height = {init_h_high}, Rise Time = {init_tr_high}, Tau = {init_tau_high}")
print(f"Initial Simulation Low Pulse : Height = {init_h_low}, Rise Time = {init_tr_low}, Tau = {init_tau_low}")

plt.plot(time_axis, pulse_high, label="Experimental High Pulse", color="gray", linestyle="--")
plt.plot(time_axis, pulse_low, label="Experimental Low Pulse", color="gray", linestyle="--")
plt.plot(time_axis, pulse_sim_high, label="Simulated High Pulse")
plt.plot(time_axis, pulse_sim_low, label="Simulated Low Pulse")
plt.xlabel("Time (s)")
plt.ylabel("Amplitude")
plt.legend()
plt.show()


progress = {
    "eval_count": 0,
    "callback_count": 0,
    "start_time": pytime.time(),
    "last_error": None,
}


def err_func(params):
    progress["eval_count"] += 1

    para["C_abs"] = params[0]
    para["C_tes"] = params[1]
    para["G_abs-abs"] = params[2]
    para["G_abs-tes"] = params[3]
    para["G_tes-bath"] = params[4]

    if any(p <= 0 for p in params):
        return 1e18

    try:
        pulses = MakePulse(para, Presample)
        pulse_sim_high = pulses[0]
        pulse_sim_low = pulses[-1]
        sim_h_low, sim_tr_low, _ = extract_features(pulse_sim_low, dt)
        sim_h_high, sim_tr_high, sim_tau = extract_features(pulse_sim_high, dt)
    except Exception:
        return 1e18

    e_wave_high = weighted_wave_error(pulse_sim_high, pulse_high, wave_weights_high, Presample)
    e_wave_low = weighted_wave_error(pulse_sim_low, pulse_low, wave_weights_low, Presample)

    e_h_high = ((sim_h_high - exp_h_high) / exp_h_high) ** 2
    e_tr_high = ((sim_tr_high - exp_tr_high) / exp_tr_high) ** 2
    e_h_low = ((sim_h_low - exp_h_low) / exp_h_low) ** 2
    e_tr_low = ((sim_tr_low - exp_tr_low) / exp_tr_low) ** 2
    e_tau = ((sim_tau - exp_tau) / exp_tau) ** 2

    wave_errors = np.array([e_wave_high, e_wave_low])
    wave_error_weights = np.array([0.0, 0.0])

    feature_errors = np.array([e_h_high, e_tr_high, e_h_low, e_tr_low, e_tau])
    feature_error_weights = np.array([8.0, 12.0, 8.0, 6.0, 6.0])

    total_error = np.sum(wave_error_weights * wave_errors) + np.sum(feature_error_weights * feature_errors)
    progress["last_error"] = total_error

    if progress["eval_count"] == 1 or progress["eval_count"] % 20 == 0:
        elapsed = pytime.time() - progress["start_time"]
        print(
            f"[Optimize] eval={progress['eval_count']}, "
            f"error={total_error:.6e}, "
            f"wave_high={e_wave_high:.3e}, wave_low={e_wave_low:.3e}, "
            f"elapsed={elapsed:.1f}s"
        )

    return total_error


def optimization_callback(xk):
    progress["callback_count"] += 1
    elapsed = pytime.time() - progress["start_time"]
    params_str = ", ".join(f"{v:.3e}" for v in xk)
    error_str = f"{progress['last_error']:.6e}" if progress["last_error"] is not None else "N/A"
    print(
        f"[Optimize] iter={progress['callback_count']}, "
        f"last_error={error_str}, elapsed={elapsed:.1f}s"
    )
    print(f"           params=[{params_str}]")


initial_params = [7.9e-10, 7.9e-12, 1.5e-07, 8.2e-09, 1.68e-08]

print("[Optimize] Starting Nelder-Mead optimization")
result = minimize(
    err_func,
    initial_params,
    method="Nelder-Mead",
    callback=optimization_callback,
    options={"maxiter": 500, "disp": True, "xatol": 1e-8, "fatol": 1e-8},
)

print("Optimization Result:")
print(result.x)

para["C_abs"], para["C_tes"], para["G_abs-abs"], para["G_abs-tes"], para["G_tes-bath"] = result.x
with open(f"{sim_path}/final_params.json", "w") as f:
    json.dump(para, f, indent=4)

pulses = MakePulse(para, Presample)
final_high = pulses[0]
final_low = pulses[-1]

final_h_high, final_tr_high, final_tau_high = extract_features(final_high, dt)
final_h_low, final_tr_low, final_tau_low = extract_features(final_low, dt)
print(f"High Pulse Simulated: Height = {final_h_high}, Rise Time = {final_tr_high}, Tau = {final_tau_high}")
print(f"Low Pulse Simulated : Height = {final_h_low}, Rise Time = {final_tr_low}, Tau = {final_tau_low}")

plt.plot(time_axis, pulse_high, label="Experimental High Pulse", color="gray", linestyle="--")
plt.plot(time_axis, pulse_low, label="Experimental Low Pulse", color="gray", linestyle="--")
plt.plot(time_axis, final_high, label="Simulated High Pulse")
plt.plot(time_axis, final_low, label="Simulated Low Pulse")
plt.xlabel("Time (s)")
plt.ylabel("Amplitude")
plt.legend()
plt.show()
