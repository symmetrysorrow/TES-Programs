import json
import matplotlib.pyplot as plt
import numpy as np
import sys

batch_path=input("batch.json:")

# JSONデータの読み込み
with open(batch_path, "r") as f:
    data = json.load(f)

event_number=input("Event number:")

if event_number in data:
    data = data[event_number]
else:
    print(f"{event_number} はデータに存在しません")
    sys.exit()

# 粒子タイプごとの色とラベル
particle_types = {
    12: ("blue", "Electron"),
    13: ("red", "Positron"),
    14: ("green", "Photon")
}

# 2Dプロットの作成
fig, ax = plt.subplots()

# 背景を鉛色にし、端を黒で囲む
ax.add_patch(plt.Rectangle((-1, -0.1), 2, 0.1, color='gray', alpha=0.5, edgecolor='black', linewidth=2))

# x軸を1~-1の範囲で300分割し、y=0~-0.1の範囲で線を引く
x_grid = np.linspace(-1, 1, 101)
for x in x_grid:
    ax.plot([x, x], [0, -0.1], color='black', linewidth=0.5, alpha=0.3)

previous_end = None

total_energy=0.0

# 軌跡をプロットし、各要素の最後の点と次の要素の最初の点をつなぐ
for key in sorted(data.keys(), key=int):
    values = data[key]
    x = values["x"]
    z = values["z"]
    ityp = values.get("ityp", 14)  # デフォルトをPhotonに設定
    color, label = particle_types.get(ityp, ("black", "Unknown"))
    
    ax.plot(x, z, marker='o', color=color, label=f'{label} ({key})')
    
    # 接続線をプロット
    if previous_end is not None:
        ax.plot([previous_end[0], x[0]], [previous_end[1], z[0]], marker='o', linestyle='--', color='gray')
    
    # 現在の要素の最後の点を保存
    previous_end = (x[-1], z[-1])

    E_deposit=sum(values["E_deposit"])
    total_energy+=E_deposit


print(f"Total deposited energy: {total_energy} MeV")

# 軸ラベル
#ax.set_xlabel("X")
#ax.set_ylabel("Z")
#ax.set_title("2D Trajectories (Y=0)")
plt.xlim(-1,1)
plt.ylim(-0.11,0.1)
ax.legend()

# 表示
plt.show()
