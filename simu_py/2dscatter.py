import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from InquirerPy import inquirer
from matplotlib.widgets import RectangleSelector

folder = input("target folder: ")

df0 = pd.read_csv(folder + "/output_TES0.csv", header=0)
df1 = pd.read_csv(folder + "/output_TES1.csv", header=0)

common_ids = set(df0['id']) & set(df1['id'])

df0 = df0[df0['id'].isin(common_ids)].reset_index(drop=True)
df1 = df1[df1['id'].isin(common_ids)].reset_index(drop=True)

columns = ['height', 'rise', 'ST_Height']
channels = [0, 1]

x_axis_ch = inquirer.select(
    message="X axis Channel:",
    choices=channels
).execute()
x_axis = inquirer.select(
    message="X axis Target:",
    choices=columns
).execute()

y_axis_ch = inquirer.select(
    message="Y axis Channel:",
    choices=channels
).execute()
y_axis = inquirer.select(
    message="Y axis Target:",
    choices=columns
).execute()

x_df = df0 if x_axis_ch == 0 else df1
y_df = df0 if y_axis_ch == 0 else df1

x_data = x_df[x_axis].values
y_data = y_df[y_axis].values

if "id" in x_df.columns:
    x_ids = x_df["id"].values
else:
    x_ids = x_df.iloc[:, 0].values

# ===== ここだけ追加 =====
valid_mask = np.isfinite(x_data) & np.isfinite(y_data)

x_data = x_data[valid_mask]
y_data = y_data[valid_mask]
x_ids  = x_ids[valid_mask]
# ======================

fig, ax = plt.subplots()
sc = ax.scatter(x_data, y_data)

selected_ids = set()

def onselect(eclick, erelease):
    global selected_ids

    x_min, x_max = sorted([eclick.xdata, erelease.xdata])
    y_min, y_max = sorted([eclick.ydata, erelease.ydata])

    newly_selected = set()

    for i in range(len(x_data)):
        if x_min <= x_data[i] <= x_max and y_min <= y_data[i] <= y_max:
            newly_selected.add(x_ids[i])

    selected_ids.update(newly_selected)

    with open(f"{folder}/selected_ids.txt", "w") as f:
        for id_val in sorted(selected_ids):
            f.write(f"{id_val}\n")

    print("Saved to selected_ids.txt")

    mask = np.isin(x_ids, list(newly_selected))
    ax.scatter(x_data[mask], y_data[mask], color='red', label='Selected')

    fig.canvas.draw()

toggle_selector = RectangleSelector(
    ax, onselect, useblit=True,
    button=[1],
    minspanx=5, minspany=5, spancoords='pixels',
    interactive=True
)

print(f"min:{np.min(y_data)},max:{np.max(y_data)}")
plt.title("2D Scatter")
plt.ylim(np.min(y_data)*0.9, np.max(y_data)*1.1)
plt.show()
