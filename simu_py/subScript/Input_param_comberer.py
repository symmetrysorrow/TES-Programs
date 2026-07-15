import json

path=input("input path:")
with open(path, "r") as f:
    input_para = json.load(f)

pre_T_c= input_para["T_c"]
pre_T_bath= input_para["T_bath"]

print(f"pre_T_c: {pre_T_c}, pre_T_bath: {pre_T_bath}")
new_T_c = float(input("input new T_c: "))
new_T_bath = float(input("input new T_bath: "))
input_para["T_c"] = new_T_c
input_para["T_bath"] = new_T_bath

new_C_abs=input_para["C_abs"]*(new_T_c/pre_T_c)**3
input_para["C_abs"] = new_C_abs

new_C_tes=input_para["C_tes"]*(new_T_c/pre_T_c)**3
input_para["C_tes"] = new_C_tes

new_G_abs_abs=input_para["G_abs-abs"]*(new_T_c/pre_T_c)**3
input_para["G_abs-abs"] = new_G_abs_abs

new_G_abs_tes=input_para["G_abs-tes"]*(new_T_c/pre_T_c)**3
input_para["G_abs-tes"] = new_G_abs_tes

new_G_tes_bath=input_para["G_tes-bath"]*(new_T_c/pre_T_c)**3
input_para["G_tes-bath"] = new_G_tes_bath

with open(path, "w") as f:
    json.dump(input_para, f, indent=4)
