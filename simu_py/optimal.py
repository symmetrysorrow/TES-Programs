import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import simpson
import sys
import tqdm
import json
import glob
from scipy.optimize import curve_fit
from scipy import signal

def GetFFT(data,rate,sample):
    g_f=np.fft.fft(data)/fs
    frequencies=np.fft.fftfreq(int(sample), d=1/rate)
    g_f_shifted = np.fft.fftshift(g_f)
    frequencies_shifted = np.fft.fftshift(frequencies)
    return g_f_shifted,frequencies_shifted
def GetIFFT(data,rate):
    g_ifft=np.fft.ifft(data)*rate
    return g_ifft

def optimal_bin_count(data):
    q1, q3 = np.percentile(data, [25, 75])
    iqr = q3 - q1  # 四分位範囲
    bin_width = 2 * iqr / (len(data) ** (1/3))  # ビン幅
    bin_count = int(np.ceil((np.max(data) - np.min(data)) / bin_width))  # ビン数
    return max(bin_count, 1)  # ビン数が1未満にならないようにする

def gaussian(x, amp, mean, stddev):
    return amp * np.exp(-((x - mean) ** 2) / (2 * stddev ** 2))

def BesselFilter(x, rate, fs):
	ws = fs / rate * 2
	b, a = signal.bessel(2, ws, "low")
	y = signal.filtfilt(b, a, x)
	return y

def MakeHistgram(data,posi,HistColor=None):
    bin_num = optimal_bin_count(data)
    bin_num=40
    hist, bin_edges = np.histogram(data, bins=bin_num, density=False)
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2  # ビンの中心を計算
    initial_guess = [np.max(hist), np.mean(data), np.std(data)]
    if HistColor is not None:
        plt.hist(data, bins=bin_num, density=False, label=f"abs-{posi}",color=HistColor)
    else:
        plt.hist(data, bins=bin_num, density=False, label=f"abs-{posi}")
    # ガウスフィッティング
    popt, _pcov = curve_fit(gaussian, bin_centers, hist, p0=initial_guess, maxfev=10000000)
    _amp_fit, mean_fit, stddev_fit = popt
    fwhm = 2 * stddev_fit * np.sqrt(2 * np.log(2))

    #ヒストグラム
    x_fit = np.linspace(bin_edges[0], bin_edges[-1], 1000)  # フィッティング用のx
    plt.plot(x_fit, gaussian(x_fit, *popt),color="red",alpha=0.5)  # フィッティング曲線

    #plt.axvline(x=mean_fit)
    #plt.show()
    return fwhm,fwhm/mean_fit

with open("f:/hata/1332_adaptive/input.json") as f:
    data = json.load(f)

for posi in data["position"]:
    for ch in [0,1]:

        target=f"f:/hata/1332_adaptive/1332keV_{posi}"

        fs = 1e6  # サンプリング周波数 (Hz)
        N=1e6
        t = np.linspace(0, 1, int(fs), endpoint=False)  # 時間軸
        # ----データの読み込み----
        M_t=np.loadtxt(f"{target}/Pulse/CH{ch}/CH{ch}_1.dat")
        N_t=np.loadtxt("f:/hata/1332_adaptive/noise_time_domain.dat")
        # ----データの読み込み----
        padded_data = np.pad(M_t, (10000, 0), mode='constant')  # 先頭に0を追加
        M_t = padded_data[:-10000]  # 最後のn個を削除
        D_t=M_t+N_t
        # ----フーリエ変換----
        M_f,frequencies = GetFFT(M_t,fs,N)
        N_f,frequencies=GetFFT(N_t,fs,N)
        D_f,frequencies=GetFFT(D_t,fs,N)
        # ----フーリエ変換----
        ratio_squared = np.abs(M_f / N_f)**2

        # 時間領域でのエネルギー計算
        energy_time = np.sum(D_t**2) * (t[1] - t[0])  # 数値積分

        # 周波数領域でのエネルギー計算
        power_spectrum_M = np.abs(D_f)**2  # パワースペクトル
        energy_freq = simpson(power_spectrum_M, x=frequencies)  # 数値積分

        # 結果の表示
        #print(f"時間領域でのエネルギー: {energy_time}")
        #print(f"周波数領域でのエネルギー: {energy_freq}")

        ratio_squared = np.abs(M_f / N_f)**2

        C=simpson(ratio_squared,x=frequencies)


        # 結果の出力
        print(f"∫|M/N|^2 df = {C}")

        # 結果のプロット
        plt.figure(figsize=(8, 6))

        plt.plot(t,D_t,label="Pulse+Noise")
        #plt.plot(t,N_t,label="Noise")
        plt.plot(t, M_t, label="Model Pulse")
        #plt.title('Time Domain Signal')
        plt.xlabel('Time [s]')
        plt.ylabel('Amplitude[A]')
        plt.xlim(0,0.2)
        plt.legend()
        plt.grid()
        plt.tight_layout()
        plt.savefig(f"{target}/{ch}_time_domain.png")#,transparent=True)
        #plt.show()
        plt.cla()
        

        #plt.plot(frequencies,np.abs(D_f),label="Pulse+Noise")
        #plt.plot(frequencies, np.abs(M_f),label="Model Pulse")
        #plt.plot(frequencies,np.abs(N_f),label="Noise")
        #plt.title('Frequency Domain Signal')
        #plt.xlabel('Frequency [Hz]')
        #plt.ylabel('Amplitude')
        #plt.grid()
        #plt.legend()
        #plt.loglog()
        #plt.tight_layout()
        #plt.savefig(f"{target}/{ch}_time_domain.png",transparent=True)
        #plt.cla()
        #plt.show()

        #Template=GetIFFT(np.conj(M_f)/((np.abs(N_f))**2),N)
        Template=GetIFFT(M_f/((np.abs(N_f))**2),N)
        Template=np.abs(Template)

        plt.plot(t,Template)
        plt.xlabel("Time[s]")
        plt.xlim(0,0.2)
        plt.savefig(f"{target}/{ch}_template.png")#,transparent=True)
        plt.show()
        plt.cla()

        pulse_noise=np.loadtxt(f"{target}/Pulse_noise/CH{ch}/CH{ch}_1.dat")
        pulse_noise_bessel=BesselFilter(pulse_noise,fs,1e4)
        pulse_noise*=Template
        pulse_noise_bessel*=Template

        files=glob.glob(f"{target}/Pulse_noise/CH{ch}/*.dat")

        Amps=[]
        Amps_bessel=[]
        Amps_bessel_both=[]
        Template_bessel=BesselFilter(Template,fs,1e4)

        #for file in tqdm.tqdm(files):
        #    pulse=np.loadtxt(file)
        #
        #    Amplitude=pulse*Template
        #    Amps.append(np.sum(Amplitude))
        #
        #    pulse=BesselFilter(pulse,fs,1e4)
        #    Amplitude_bessel=pulse*Template
        #   Amps_bessel.append(np.sum(Amplitude_bessel))

        #   Amplitude_bessel_both=pulse*Template_bessel
        #   Amps_bessel_both.append(np.sum(Amplitude_bessel_both))

        #np.savetxt(f"{target}/{ch}.txt",Amps)

        Amps=np.array(Amps)
        Amps_bessel=np.array(Amps_bessel)
        Amps_bessel_both=np.array(Amps_bessel_both)

        #fwhm,reso=MakeHistgram(Amps,17)
        #fwhm_bessel,reso_bessel=MakeHistgram(Amps_bessel,18)
        #fwhm_both,reso_both=MakeHistgram(Amps_bessel_both,19)
        #plt.legend()
        plt.savefig(f"{target}/{ch}_plot.png",transparent=True)
        plt.cla()
        #plt.show()

        #print(f"reso: {reso*1332}")
        #print(f"reso_bessel: {reso_bessel*1332}")
        #print(f"reso_both: {reso_both*1332}")