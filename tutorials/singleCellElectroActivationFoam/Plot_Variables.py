import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# ref_data = np.genfromtxt('bana.txt', names=True)
# print(ref_data.dtype.names)  # Check the column names
# print(ref_data.dtype.names)  # should print ('time', 'cellv')

# # Access columns by field names
# ref_time = ref_data['time']
# ref_voltage = ref_data['cellv']

# --- Load simulation output with headers ---
file = 'banana22_output.txt'
with open(file, 'r') as f:
    header = f.readline().strip().split()

data = np.loadtxt(file, skiprows=1)
df1 = pd.DataFrame(data, columns=header)

# Used for labeling the dataset
file_list = [file]

# --- Plotting function ---
def plot_variables(df_list, file_list, var_names, y_labels, figsize=(12, 12)):
    time_list = [df["time"] for df in df_list]

    fig, axes = plt.subplots(4, 2, figsize=figsize, sharex=True)
    axes = axes.flatten()
    custom_labels = [f.replace('_output.txt', '') for f in file_list]

    for ax, var, label in zip(axes, var_names, y_labels):
        for df, t, run_label in zip(df_list, time_list, custom_labels):
            ax.plot(t, df[var], label=run_label)
        ax.set_ylabel(label)
        ax.grid(True)
        ax.legend()

    for ax in axes[-2:]:
        ax.set_xlabel("Time (ms)")

    plt.tight_layout()
    plt.show()

# --- First plot: main variables ---
var_names_1 = [
    "cell_v",
    "ionic_concentrations_cai",
    "AV_INa",
    "AV_ITo_ITo",
    "AV_IKr_IKr",
    "AV_IKs_IKs",
    "AV_ICaL_ICaL",
    "AV_IK1_IK1"
]
var_names_1 = [
    "V" 
]

y_labels_1 = [
    "Voltage (mV)",
    "Ca_i (mM)",
    "i_Na (pA/pF)",
    "i_to (pA/pF)",
    "i_Kr (pA/pF)",
    "i_Ks (pA/pF)",
    "i_CaL (pA/pF)",
    "i_K1 (pA/pF)"
]
y_labels_1 = [
    "Voltage (mV)"
]

plot_variables([df1], file_list, var_names_1, y_labels_1)

# --- Plot reference voltage vs simulation voltage ---
plt.figure(figsize=(10, 4))
plt.plot(ref_time, ref_voltage, label='Reference Voltage', linestyle='--')
plt.plot(df1["time"], df1["cell_v"], label='Simulated Voltage', alpha=0.8)
plt.xlabel("Time (ms)")
plt.ylabel("Voltage (mV)")
plt.title("Simulated vs Reference Voltage")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()
