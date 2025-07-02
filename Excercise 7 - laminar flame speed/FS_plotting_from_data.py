import pandas as pd
import numpy as np

# needed to save the figures without displaying them
#import matplotlib
#matplotlib.use('Agg')  # Use Agg backend for non-interactive plotting

import matplotlib.pyplot as plt
import cantera as ct

# -------------------- Experimental Data import Section --------------------------------
# Import experimental flame speed data
exp_const_temp = pd.read_csv('/home/brukare/Combustion/Excercise 7 - laminar flame speed/FlameSpeedsExperimental_constant temperature_298K.csv')
exp_const_P = pd.read_csv('/home/brukare/Combustion/Excercise 7 - laminar flame speed/FlameSpeedsExperimental_constant pressure_1atm.csv')

# Only take coloumns with relevant data
exp_const_temp = exp_const_temp.iloc[:, [0, 1, 2, 3, 8, 9, 12, 13]]
exp_const_P = exp_const_P.iloc[:, :-3]


# -------------------- Cantera Calculation import Section  --------------------------------
# Simulation Parameters
pressure = [101325, 202650, 506635, 2026500]  # pressure in Pa]
Tin_const = 298.0  # unburned gas temperature in K
Tin_var = [358.0, 393.0, 428.0]  # unburned gas temperature in K
width = 0.03  # width of the flame in m

phi = np.linspace(0.5, 1.5, 21)  # equivalence ratio range
loglevel = 1  # amount of diagnostic outptut (0 to 8)

fin_P_Result=pd.read_csv('/home/brukare/Combustion/Excercise 7 - laminar flame speed/fin_P_Result.csv')
fin_Temp_Result=pd.read_csv('/home/brukare/Combustion/Excercise 7 - laminar flame speed/fin_Temp_Result.csv')



# -------------------- Plot Section --------------------------------

# Marker styles and colors to differentiate datasets
markers = ['o', 's', '^', 'D', 'v', '*', 'x', 'P']
colors = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple',
          'tab:brown', 'tab:gray']


fig, axs = plt.subplots(2, 1, figsize=(10, 12))

def data_to_plot(exp_data, markers, colors, ax, labels):
    num_columns = exp_data.shape[1]
    for i in range(0, num_columns, 2):
        x = exp_data.iloc[:, i]
        y = exp_data.iloc[:, i + 1]
        idx = i // 2
        marker = markers[idx % len(markers)]
        color = colors[idx % len(colors)]
        label = labels[idx % len(labels)] + " (Exp)" 
        ax.plot(x, y, linestyle='None', marker=marker, color=color, label=label)


# Plot first dataset (top subplot)
labels_temp = ["358 k", "393 k", "428 k"]
data_to_plot(exp_const_P, markers, colors, axs[0], labels_temp)
for i in range(len(Tin_var)):
    plots = axs[0].plot(phi, fin_Temp_Result.iloc[i, 1:], linestyle='-', color=colors[i], label= labels_temp[i] + ' (Cantera)')
axs[0].set_title('Flame Speed at Constant Pressure (1 atm)', fontsize=18)
axs[0].set_xlabel('Equivalence Ratio ϕ', fontsize=16)
axs[0].set_ylabel('Laminar Flame Speed [cm/s]', fontsize=16)
axs[0].legend(fontsize=14)
axs[0].tick_params(axis='both', labelsize=14)
axs[0].grid(True)

# Plot second dataset (bottom subplot)
labels_presure = ["1 atm", "2 atm", "5 atm", "20 atm"]
data_to_plot(exp_const_temp, markers, colors, axs[1], labels_presure)
for i in range(len(pressure)):
    plots = axs[1].plot(phi, fin_P_Result.iloc[i, 1:], linestyle='-', color=colors[i], label= labels_presure[i] + ' (Cantera)')
axs[1].set_title('Flame Speed at Constant Temperature (298 K)', fontsize=18)
axs[1].set_xlabel('Equivalence Ratio ϕ', fontsize=16)
axs[1].set_ylabel('Laminar Flame Speed [cm/s]', fontsize=16)
axs[1].legend(fontsize=14)
axs[1].tick_params(axis='both', labelsize=14)
axs[1].grid(True)

plt.tight_layout()
plt.show()

# Save the figure, won't display it, need to activate Agg backend under imports
#plt.savefig('/home/brukare/Combustion/Excercise 7 - laminar flame speed/FlameSpeed_Results_Both.png', dpi=300)
