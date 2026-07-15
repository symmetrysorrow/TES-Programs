import os
import json

Output = "F:/hata/1332_142_136_500split"

with open(Output + "/input.json", 'r') as file:
    data = json.load(file)

positions = data['position']
output_keV = int(data["E"])
history = data['history']  # 直接整数として扱う
n_abs=data["n_abs"]
length=data["length"]

with open(Output + "/trace.inp", 'r') as file:
    trace_content = file.readlines()

# x0= * と maxcas= * を置換する関数
def replace_trace_content(k, trace_content):
    new_content = []
    for line in trace_content:
        # x0= * の行を書き換える
        if line.startswith('x0'):
            x0_value = length*(2*k-n_abs-1)/(2*n_abs*10)
            new_x0 = f"x0 = {x0_value}"
            new_content.append(new_x0)
        # maxcas= * の行をhistoryに書き換える
        elif line.startswith('maxcas'):
            new_maxcas = f"maxcas = {int(history)}"
            new_content.append(new_maxcas)
        else:
            new_content.append(line.strip())  # 改行を除去しておく
    return new_content

# 各positionに対してフォルダ作成とtrace.inpの書き換えを行う
for pos in positions:
    folder_name = f"{output_keV}keV_{pos}"
    folder_path = os.path.join(Output, folder_name)

    # フォルダがない場合は作成
    os.makedirs(folder_path, exist_ok=True)

    modified_trace_content = replace_trace_content(pos, trace_content)

    # 新しいtrace.inpをフォルダ内に書き込む
    output_trace_path = os.path.join(folder_path, 'trace.inp')
    with open(output_trace_path, 'w') as file:
        for line in modified_trace_content:
            print(line, file=file)  # print() を使って改行を自動的に挿入

print("trace.inpファイルが各フォルダに作成されました。")
